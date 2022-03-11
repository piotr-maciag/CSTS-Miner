#include "load_data.h"

double T;
double R;


//########################################

double theta;
double approxMargin;

double phase1Time;
double phase2Time;

int sizeP;
STPoint *dataP;
fstream uchwyt;

int numOfEventTypes; // not used

vector<vector<STPoint>> dataset;
vector<vector<STPoint>> sortedDataset;

vector<vector<Sequence *>> SequencesSetSPTree;

int GseqID = 0;


bool comparison(const STPoint &i1, const STPoint &i2) {
    return i1.temporal < i2.temporal;
}

bool comparisonPT(STPoint *i1, STPoint *i2) {
    return i1->temporal < i2->temporal;
}

bool comparisonID(const STPoint &i1, const STPoint &i2) {
    return i1.eventID < i2.eventID;
}

bool comparisonIDPT(STPoint *i1, STPoint *i2) {
    return i1->eventID < i2->eventID;
}

bool isEqual(const STPoint &i1, const STPoint &i2) {
    return (i1.eventID == i2.eventID);
}

bool isNotEqual(const STPoint &i1, const STPoint &i2) {
    return (i1.eventID != i2.eventID);
}

void SortDataset() {
    sortedDataset = dataset;
    for (int i = 0; i < dataset.size(); i++) {
        sort(sortedDataset[i].begin(), sortedDataset[i].end(), comparison);
    }
}


#define earthRadiusKm 6371.0

// This function converts decimal degrees to radians
double deg2rad(double deg) {
    return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
    return (rad * 180 / M_PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
    double lat1r, lon1r, lat2r, lon2r, u, v;
    lat1r = deg2rad(lat1d);
    lon1r = deg2rad(lon1d);
    lat2r = deg2rad(lat2d);
    lon2r = deg2rad(lon2d);
    u = sin((lat2r - lat1r) / 2);
    v = sin((lon2r - lon1r) / 2);
    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v)) * 1000;
}

double distanceSpatial(double pSpatialY, double pSpatialX, double qSpatialY, double qSpatialX) {
    return sqrt((pSpatialX - qSpatialX) * (pSpatialX - qSpatialX) + (pSpatialY - qSpatialY) * (pSpatialY - qSpatialY));
}

vector<STPoint *> ForwardSweep(vector<STPoint *> tailEventSet, vector<STPoint *> instancesSet) {

    sort((tailEventSet.begin()), (tailEventSet.end()), comparisonPT);
    sort(instancesSet.begin(), instancesSet.end(), comparisonPT);
    vector<STPoint *> joinResult;


    while (tailEventSet.empty() != true && instancesSet.empty() != true) {
        int pindex = 0, qindex = 0;

        STPoint *p = tailEventSet[pindex];
        STPoint *q = instancesSet[qindex];

        if (p->temporal < q->temporal) {
            tailEventSet.erase(tailEventSet.begin());
            while (p->temporal + T > q->temporal) {

                double dist = distanceEarth(p->spatialY, p->spatialX, q->spatialY, q->spatialX);

                if (dist <= R) {
                    joinResult.push_back(q);
                }
                if (qindex < instancesSet.size() - 1) {
                    qindex++;
                    q = instancesSet[qindex];
                } else {
                    break;
                }
            }
        } else {
            instancesSet.erase(instancesSet.begin());
        }
    }

    std::vector<int>::iterator it;
//	sort( joinResult.begin(), joinResult.end(), comparisonID);
    //joinResult.erase(unique(joinResult.begin(), joinResult.end(), isEqual), joinResult.end());

    sort(joinResult.begin(), joinResult.end(), comparisonIDPT);
    joinResult.erase(unique(joinResult.begin(), joinResult.end()), joinResult.end());


    return joinResult;
}

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

double CalculatePI(vector<STPoint *> instSet) {
    double totalNumber;

    if (instSet.empty() == true) {
        return 0;
    }

    for (int i = 0; i < sortedDataset.size(); i++) {
        if (instSet[0]->eventType == sortedDataset[i][0].eventType) {
            totalNumber = sortedDataset[i].size();
        }
    }


    return (double(instSet.size()) / totalNumber);
}



//##################################################################
//##################################################################

//void RemoveFromRC(Sequence *s) {
//    for (int i = 0; i < s->Cmax->RC.size(); i++) {
//        if (s->Cmax->RC[i] == s) {
//            s->Cmax->RC.erase(s->Cmax->RC.begin() + i);
//        }
//    }
//}

void CSTS_Miner() {

    vector<Sequence *> Length1Seq;
    vector<Sequence *> Length2Seq;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (int i = 0; i < sortedDataset.size(); i++) //create 1-length sequences
    {
        Sequence *seq = new Sequence;
        seq->sequence.push_back(sortedDataset[i][0].eventType);
        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }
        seq->tailEventSet.push_back(pointersData);
        GseqID++;
        //seq->seqID = GseqID;
        seq->PI = 1.0;
        seq->Cmax.clear();
        seq->RC.clear();
        Length1Seq.push_back(seq);
    }
    SequencesSetSPTree.push_back(Length1Seq);
    cout << "#" << flush;

    for (int i = 0; i < Length1Seq.size(); i++) //create 2-length sequences
    {
        for (int j = 0; j < Length1Seq.size(); j++) {
            Sequence *seq = new Sequence;

            seq->sequence.push_back(Length1Seq[i]->sequence[0]);
            seq->sequence.push_back(Length1Seq[j]->sequence[0]);
            vector<STPoint *> I2 = ForwardSweep(Length1Seq[i]->tailEventSet[0], Length1Seq[j]->tailEventSet[0]);

            if (I2.empty() == false) {
                seq->tailEventSet.push_back(I2);
                seq->PI = CalculatePI(seq->tailEventSet[0]);

                if (seq->PI > theta) {
                    GseqID++;
                    //seq->seqID = GseqID;
                    seq->firstParent = Length1Seq[i];
                    seq->secondParent = Length1Seq[j];
                    seq->firstParent->children.push_back(seq);
                    seq->Cmax.clear();
                    seq->RC.clear();
                    Length2Seq.push_back(seq);
                } else
                {
                    delete seq;
                }
            }
        }
    }

    SequencesSetSPTree.push_back(Length2Seq);
    cout << "#" << flush;

    int k = 1;
    while (SequencesSetSPTree[k].empty() == false) {
        k++;
        vector<Sequence *> LengthKSeq = GenAndVerify(SequencesSetSPTree[k - 1]);
        SequencesSetSPTree.push_back(LengthKSeq);
        cout << "#" << flush;
        for(int i = 0; i < SequencesSetSPTree[k-1].size(); i++)
        {
            SequencesSetSPTree[k-1][i]->tailEventSet.clear();
        }
    }


    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double duration1 = duration_cast<milliseconds>(t2 - t1).count();

    cout << endl;

    t1 = high_resolution_clock::now();
    while (k > 0) {
        for (int i = 0; i < SequencesSetSPTree[k].size(); i++) {
            VerifySupersequence(SequencesSetSPTree[k][i], SequencesSetSPTree[k][i]->firstParent);
            VerifySupersequence(SequencesSetSPTree[k][i], SequencesSetSPTree[k][i]->secondParent);
        }
        cout << "&"<< flush;
        k--;
    }
    t2 = high_resolution_clock::now();
    double duration2 = duration_cast<milliseconds>(t2 - t1).count();

    phase1Time = duration1/1000;
    phase2Time = duration2/1000;

    cout << endl;
}

void STBFM() {

    vector<Sequence *> Length1Seq;
    vector<Sequence *> Length2Seq;
    GseqID = 0;

    for (int i = 0; i < sortedDataset.size(); i++) //create 1-length sequences
    {
        Sequence *seq = new Sequence;
        seq->sequence.push_back(sortedDataset[i][0].eventType);
        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }
        seq->tailEventSet.push_back(pointersData);
        GseqID++;
        //seq->seqID = GseqID;
        seq->PI = 1.0;
        seq->Cmax.clear();
        seq->RC.clear();
        Length1Seq.push_back(seq);
    }
    SequencesSetSPTree.push_back(Length1Seq);

    for (int i = 0; i < Length1Seq.size(); i++) //create 2-length sequences
    {
        for (int j = 0; j < Length1Seq.size(); j++) {
            Sequence *seq = new Sequence;

            seq->sequence.push_back(Length1Seq[i]->sequence[0]);
            seq->sequence.push_back(Length1Seq[j]->sequence[0]);
            vector<STPoint *> I2 = ForwardSweep(Length1Seq[i]->tailEventSet[0], Length1Seq[j]->tailEventSet[0]);

            if (I2.empty() == false) {
                seq->tailEventSet.push_back(I2);
                seq->PI = CalculatePI(seq->tailEventSet[0]);

                if (seq->PI > theta) {
                    GseqID++;
                    //seq->seqID = GseqID;
                    seq->firstParent = Length1Seq[i];
                    seq->secondParent = Length1Seq[j];
                    seq->firstParent->children.push_back(seq);
                    Length2Seq.push_back(seq);
                }else
                {
                    delete seq;
                }
            }
        }
    }

    SequencesSetSPTree.push_back(Length2Seq);


    int k = 1;
    while (SequencesSetSPTree[k].empty() == false) {
        k++;
        vector<Sequence *> LengthKSeq = GenAndVerify(SequencesSetSPTree[k - 1]);
        SequencesSetSPTree.push_back(LengthKSeq);
        for(int i = 0; i < SequencesSetSPTree[k-1].size(); i++)
        {
            SequencesSetSPTree[k-1][i]->tailEventSet.clear();
        }
        //cout << SequencesSetSPTree[k].size() << endl;
    }
}

void CST_SPMiner() {

    vector<Sequence *> Length1Seq;
    vector<Sequence *> Length2Seq;
    GseqID = 0;

    for (int i = 0; i < sortedDataset.size(); i++) //create 1-length sequences
    {
        Sequence *seq = new Sequence;
        seq->sequence.push_back(sortedDataset[i][0].eventType);
        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }
        seq->tailEventSet.push_back(pointersData);
        GseqID++;
       // seq->seqID = GseqID;
        seq->PI = 1.0;
        seq->Cmax.clear();
        seq->RC.clear();
        Length1Seq.push_back(seq);
    }
    SequencesSetSPTree.push_back(Length1Seq);
    cout << "#" << flush;

    for (int i = 0; i < Length1Seq.size(); i++) //create 2-length sequences
    {
        for (int j = 0; j < Length1Seq.size(); j++) {
            Sequence *seq = new Sequence;

            seq->sequence.push_back(Length1Seq[i]->sequence[0]);
            seq->sequence.push_back(Length1Seq[j]->sequence[0]);
            vector<STPoint *> I2 = ForwardSweep(Length1Seq[i]->tailEventSet[0], Length1Seq[j]->tailEventSet[0]);

            if (I2.empty() == false) {
                seq->tailEventSet.push_back(I2);
                seq->PI = CalculatePI(seq->tailEventSet[0]);

                if (seq->PI > theta) {
                    GseqID++;
                //    seq->seqID = GseqID;
                    seq->firstParent = Length1Seq[i];
                    seq->secondParent = Length1Seq[j];
                    seq->firstParent->children.push_back(seq);
                    seq->Cmax.clear();
                    seq->RC.clear();
                    Length2Seq.push_back(seq);
                }else
                {
                    delete seq;
                }
            }
        }
    }

    SequencesSetSPTree.push_back(Length2Seq);
    cout << "#" << flush;


    int k = 1;
    while (SequencesSetSPTree[k].empty() == false) {
        k++;
        vector<Sequence *> LengthKSeq = GenAndVerify(SequencesSetSPTree[k - 1]);
        SequencesSetSPTree.push_back(LengthKSeq);
        cout << "#" << flush;
        for(int i = 0; i < SequencesSetSPTree[k-1].size(); i++)
        {
            SequencesSetSPTree[k-1][i]->tailEventSet.clear();
        }
        //cout << SequencesSetSPTree[k].size() << endl;
    }

    cout << endl;

}

void STS_Miner()
{
    vector<Sequence *> Length1Seq;
    GseqID = 0;

    for (int i = 0; i < sortedDataset.size(); i++) //create 1-length sequences
    {
        Sequence *seq = new Sequence;
        seq->sequence.push_back(sortedDataset[i][0].eventType);
        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }
        seq->tailEventSet.push_back(pointersData);
        GseqID++;
        // seq->seqID = GseqID;
        seq->PI = 1.0;
        seq->Cmax.clear();
        seq->RC.clear();
        Length1Seq.push_back(seq);
    }
    SequencesSetSPTree.push_back(Length1Seq);
    cout << "#"<< flush;

    for(int i = 0; i < Length1Seq.size(); i++)
    {
        ExpandSequence(Length1Seq[i]);
        cout << "#" << flush;
    }
}

void ExpandSequence(Sequence * seq)
{
    int k = seq->sequence.size();

    //cout << k << " " << endl;
    int actualSize = SequencesSetSPTree.size();
    vector<Sequence *> LengthKSeq;

    if(k+1 != actualSize)
    {
        SequencesSetSPTree.push_back(LengthKSeq);
    }

    for (int i = 0; i < sortedDataset.size(); i++)
    {
        Sequence *seqNew = new Sequence;

        for (int l = 0; l < seq->sequence.size(); l++) {
            seqNew->sequence.push_back(seq->sequence[l]);
        }
        seqNew->sequence.push_back(sortedDataset[i][0].eventType);

        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }

        seqNew->tailEventSet.push_back(ForwardSweep(seq->tailEventSet[0], pointersData));

        double PI = CalculatePI(seqNew->tailEventSet[0]);
        seqNew->PI = PI;



        if (seqNew->PI > theta) {

            GseqID++;
            //seq->seqID = GseqID;

            SequencesSetSPTree[SequencesSetSPTree.size() - 1].push_back(seqNew);
            ExpandSequence(seqNew);

        }else
        {
            delete seqNew;
        }
    }
}

void CSTPM() {

    vector<Sequence *> Length1Seq;
    vector<Sequence *> Length2Seq;
    GseqID = 0;

    for (int i = 0; i < sortedDataset.size(); i++) //create 1-length sequences
    {
        Sequence *seq = new Sequence;
        seq->sequence.push_back(sortedDataset[i][0].eventType);
        vector<STPoint *> pointersData;
        for (int j = 0; j < sortedDataset[i].size(); j++) {
            pointersData.push_back(&sortedDataset[i][j]);
        }
        seq->tailEventSet.push_back(pointersData);
        GseqID++;
        //seq->seqID = GseqID;
        seq->PI = 1.0;
        seq->Cmax.clear();
        seq->RC.clear();
        Length1Seq.push_back(seq);
    }
    SequencesSetSPTree.push_back(Length1Seq);
    cout << "#" << flush;

    for (int i = 0; i < Length1Seq.size(); i++) //create 2-length sequences
    {
        for (int j = 0; j < Length1Seq.size(); j++) {
            Sequence *seq = new Sequence;


            seq->sequence.push_back(Length1Seq[i]->sequence[0]);
            seq->sequence.push_back(Length1Seq[j]->sequence[0]);
            vector<STPoint *> I2 = ForwardSweep(Length1Seq[i]->tailEventSet[0], Length1Seq[j]->tailEventSet[0]);

            if (I2.empty() == false) {
                seq->tailEventSet.push_back(I2);
                seq->PI = CalculatePI(seq->tailEventSet[0]);

                if (seq->PI > theta) {
                    GseqID++;
                    //seq->seqID = GseqID;
                    Length2Seq.push_back(seq);
                }else
                {
                    delete seq;
                }
            }
        }
    }

    SequencesSetSPTree.push_back(Length2Seq);
    cout << "#" << flush;

    int k = 1;
    while (SequencesSetSPTree[k].empty() == false) {
        k++;
        vector<Sequence *> LengthKSeq = GenAndVerifyCSTPM(SequencesSetSPTree[k - 1]);
        SequencesSetSPTree.push_back(LengthKSeq);
        cout << "#" << flush;
        for(int i = 0; i < SequencesSetSPTree[k-1].size(); i++)
        {
            SequencesSetSPTree[k-1][i]->tailEventSet.clear();
        }
        //cout << SequencesSetSPTree[k].size() << endl;
    }
}



//void RemoveSeqFromRC(Sequence *seqPar) {
//    Sequence *si;
//    for (int i = 0; i < seqPar->Cmax.size(); i++) {
//        si = seqPar->Cmax[i];
//        for (int j = 0; j < si->RC.size(); j++) {
//            if (seqPar->seqID == si->RC[j]->seqID) {
//                si->RC.erase(si->RC.begin() + j);
//            }
//        }
//    }
//}

void RemoveSeqFromRC_sl(Sequence *sl) {
    Sequence *si;
    for (int i = 0; i < sl->Cmax.size(); i++) {
        si = sl->Cmax[i];
        for (int j = 0; j < si->RC.size(); j++) {
            if (sl == si->RC[j]) {
                si->RC.erase(si->RC.begin() + j);
            }
        }
    }
}

bool belongsSeq(Sequence * seqPar, Sequence * seq)
{

    for(int i = 0; i < seq->RC.size(); i++)
    {
        //auto it = search(seq->RC[i]->sequence.begin(), seq->RC[i]->sequence.end(), boyer_moore_searcher(seqPar->sequence.begin(), seqPar->sequence.end()));

        //bool isPresent = false;

        if(seqPar->sequence.size() == seq->RC[i]->sequence.size()) {
            for (int j = 0; j < seqPar->sequence.size(); j++) {
                if (seqPar->sequence[j] == seq->RC[i]->sequence[j])
                {
                    if(j == (seqPar->sequence.size() - 1))
                    {
                        return true;
                    }
                } else
                {
                    break;
                }
            }
        }

//        if(it != seq->RC[i]->sequence.end())
//        {
//            return true;
//        }
    }

    return false;
}

void VerifySupersequence(Sequence *seq, Sequence *seqPar) {

    if (seq->PI >= seqPar->PI - approxMargin) {
        if(belongsSeq(seqPar, seq) == false) {
            if (seqPar->Cmax.empty() == true || seq->sequence.size() == seqPar->Cmax[0]->sequence.size()) {
                if (seqPar->Cmax.empty() == true || seq->PI == seqPar->Cmax[0]->PI) {
                    seqPar->Cmax.push_back(seq);
                    seq->RC.push_back(seqPar);
                } else if (seq->PI > seqPar->Cmax[0]->PI) {
                    RemoveSeqFromRC_sl(seqPar);
                    seqPar->Cmax.clear();
                    seqPar->Cmax.push_back(seq);
                    seq->RC.push_back(seqPar);
                }
            }
        }

        if (seqPar->firstParent != NULL) {
            VerifySupersequence(seq, seqPar->firstParent);
        }
        if (seqPar->secondParent != NULL) {
            VerifySupersequence(seq, seqPar->secondParent);
        }
    }
}

vector<Sequence *> GenAndVerifyCSTPM(vector<Sequence *> SeqVec) {
    vector<Sequence *> newSeqs;

    for (int i = 0; i < SeqVec.size(); i++) {
        for (int j = 0; j < SeqVec.size(); j++) {
            bool toMerge = true;

            for (int l = 0; l < SeqVec[i]->sequence.size() - 1; l++) {
                if (SeqVec[i]->sequence[l + 1] != SeqVec[j]->sequence[l])
                    toMerge = false;
            }

            if(toMerge == true) {

                Sequence *seq = new Sequence;

                for (int l = 0; l < SeqVec[i]->sequence.size(); l++) {
                    seq->sequence.push_back(SeqVec[i]->sequence[l]);
                }
                seq->sequence.push_back(SeqVec[j]->sequence[SeqVec[j]->sequence.size() - 1]);

                seq->tailEventSet.push_back(ForwardSweep(SeqVec[i]->tailEventSet[0],
                                                         SeqVec[j]->tailEventSet[0]));


                double PI = CalculatePI(seq->tailEventSet[0]);
                seq->PI = min(PI, SeqVec[i]->PI);

                if (seq->PI > theta) {

                    GseqID++;
                    newSeqs.push_back(seq);

                } else {
                    delete seq;
                }
            }
        }
    }


    return newSeqs;
}

vector<Sequence *> GenAndVerify(vector<Sequence *> SeqVec) {
    vector<Sequence *> newSeqs;

    for (int i = 0; i < SeqVec.size(); i++) {
        Sequence *sl = SeqVec[i]->secondParent;
        for (int j = 0; j < sl->children.size(); j++) {
            Sequence *seq = new Sequence;

            for (int l = 0; l < SeqVec[i]->sequence.size(); l++) {
                seq->sequence.push_back(SeqVec[i]->sequence[l]);
            }
            seq->sequence.push_back(sl->children[j]->sequence[sl->children[j]->sequence.size() - 1]);

            seq->tailEventSet.push_back(ForwardSweep(SeqVec[i]->tailEventSet[0],
                                                     sl->children[j]->tailEventSet[0]));

            seq->firstParent = SeqVec[i];
            seq->secondParent = sl->children[j];

            double PI = CalculatePI(seq->tailEventSet[0]);
            seq->PI = min(PI, SeqVec[i]->PI);

            if (seq->PI > theta) {

                GseqID++;
                //seq->seqID = GseqID;
                seq->firstParent->children.push_back(seq);

                if (seq->PI == seq->firstParent->PI) seq->firstParent->isClosed = false;
                if (seq->PI == seq->secondParent->PI) seq->secondParent->isClosed = false;

                newSeqs.push_back(seq);

            }else
            {
                delete seq;
            }

        }
    }


    return newSeqs;
}

//##################################################################
//##################################################################

void LoadDataset(string Path) //za�aduj zbi�r danych
{
    sizeP = CountInstances(Path); // zlicz l. instancji w pliku

    dataP = new STPoint[sizeP]; // rozszerz data

    uchwyt.open(Path);
    for (int i = 0; i < sizeP; i++) {
        string line;
        getline(uchwyt, line);
        stringstream linestream(line);
        string dataPortion;

        if (line != "") {
            getline(linestream, dataPortion, ' ');
            dataP[i].eventID = stoi(dataPortion);
            getline(linestream, dataPortion, ' ');
            dataP[i].eventType = dataPortion;
            getline(linestream, dataPortion, ' ');
            dataP[i].spatialX = stod(dataPortion);
            getline(linestream, dataPortion, ' ');
            dataP[i].spatialY = stod(dataPortion);
            getline(linestream, dataPortion, ' ');
            dataP[i].temporal = stod(dataPortion);
        }

        //cout << i << endl;
    }
    uchwyt.close();
}

void TransformData() {
    for (int i = 0; i < sizeP; i++) {
        InsertInstance(dataP[i]);
    }

    //delete dataP;
}

void InsertInstance(STPoint instance) {
    for (int i = 0; i < dataset.size(); i++) {
        if (dataset[i][0].eventType == instance.eventType) {
            dataset[i].push_back(instance);
            return;
        }
    }

    vector<STPoint> vect;
    vect.push_back(instance);
    dataset.push_back(vect);
}

void PrintDataset() {

    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < dataset[i].size(); j++) {
            cout << dataset[i][j].eventID << ' ' << dataset[i][j].eventType << '\t';
            //cout << i;
        }

        cout << endl;
    }
}

void PrintSortedDataset() {

    for (int i = 0; i < dataset.size(); i++) {
        for (int j = 0; j < dataset[i].size(); j++) {
            cout << sortedDataset[i][j].temporal << ' ' << sortedDataset[i][j].eventType << '\t';
            //cout << i;
        }

        cout << endl;
    }

}

void PrintSequences() {
    fstream results;
    results.open("..//ExamplesEventSeq.txt", fstream::out);

    for (int i = SequencesSetSPTree.size() - 1; i >= 0; i--) {
        for (int j = 0; j < SequencesSetSPTree[i].size(); j++) {
            Sequence *seq = SequencesSetSPTree[i][j];

            if (seq->RC.empty() == false || (seq->Cmax.empty() == true && seq->RC.empty() == true)) {
                int count = 0;
                double aPI = seq->PI;
                //cout << "i " << seq->PI << " ";
                results << i << " " << seq << " " << seq->PI << " ";


                for (int i = 0; i < seq->sequence.size(); i++) {
                    results << seq->sequence[i] << "->";
                }

                results << seq->RC.size();

                results << " || ";

                for (int l = 0; l < seq->RC.size(); l++) {
                    Sequence *src = seq->RC[l];
                    for (int k = 0; k < src->sequence.size(); k++) {
                        results << src->sequence[k] << "->";
                    }
                    results << " ::: ";
                }

                /*
                while (seq != NULL) {
                    if (seq->PI == aPI) {
                        count++;
                    }
                    //cout << seq->sequence[0] << " <- " ;
                    results << seq->sequence[0] << " <- ";
                    seq = seq->firstParent;
                }*/

                //cout << "Closed" << endl;
                results << endl;
            }
        }
    }
    results.close();
}

void ClearStructures() {

    for(int i = 0; i < dataset.size(); i++)
    {
        dataset[i].clear();
    }
    dataset.clear();
    for(int i = 0; i < sortedDataset.size(); i++)
    {
        sortedDataset[i].clear();
    }
    sortedDataset.clear();
    //delete dataP;

//    for(int i = 0; i < sizeP; i++)
//    {
//        delete (dataP +i);
//    }

    for (int i = 0; i < SequencesSetSPTree.size(); i++) {
        for (int j = 0; j < SequencesSetSPTree[i].size(); j++) {
            SequencesSetSPTree[i][j]->RC.clear();
            SequencesSetSPTree[i][j]->Cmax.clear();
            SequencesSetSPTree[i][j]->tailEventSet.clear();
            SequencesSetSPTree[i][j]->sequence.clear();
            delete SequencesSetSPTree[i][j];
        }
        SequencesSetSPTree[i].clear();
    }
    SequencesSetSPTree.clear();
}

int CountInstances(string path) // zlicz liczb� instancji w pliku
{
    uchwyt.open(path);
    string line;

    int numInstances = 0;

    while (uchwyt.eof() != true) {
        getline(uchwyt, line);

        if (line != "") {
            numInstances++;
        }
    }

    uchwyt.close();
    return numInstances;
}

int CountPatternsCSTS_Miner() {
    int count = 0;
    for (int i = 1; i < SequencesSetSPTree.size(); i++) {
        for (int j = 0; j < SequencesSetSPTree[i].size(); j++) {
            Sequence *seq = SequencesSetSPTree[i][j];
            if (seq->RC.empty() == false || (seq->Cmax.empty() == true && seq->RC.empty() == true)) {
                count++;
            }
        }
    }

    return count;
}

int CountPatternsClosed() {
    int count = 0;
    for (int i = 1; i < SequencesSetSPTree.size(); i++) {
        for (int j = 0; j < SequencesSetSPTree[i].size(); j++) {
            Sequence *seq = SequencesSetSPTree[i][j];
            if (seq->isClosed == true) {
                count++;
            }
        }
    }

    return count;
}

int CountPatterns() {
    int count = 0;
    for (int i = 1; i < SequencesSetSPTree.size(); i++) {
        for (int j = 0; j < SequencesSetSPTree[i].size(); j++) {
            count++;
        }
    }

    return count;
}

void WriteStatistics(fstream &fileRCsize, fstream &fileCmaxsize)
{

    int maxRC = 0;

    for(int i = 0; i < SequencesSetSPTree.size() - 1; i++)
    {
        if(maxRC < SequencesSetSPTree[i].size())
        {
            maxRC = SequencesSetSPTree[i].size();
        }
    }


    for(int i = 0; i < maxRC - 1; i++)
    {
        for(int j = 0; j < SequencesSetSPTree.size() - 1; j++)
        {

            if(i < SequencesSetSPTree[j].size())
            {
                Sequence * seq = SequencesSetSPTree[j][i];
                if (seq->RC.empty() == false || (seq->Cmax.empty() == true && seq->RC.empty() == true)) {
                    fileRCsize << SequencesSetSPTree[j][i]->RC.size() << ",";
                    fileCmaxsize << SequencesSetSPTree[j][i]->Cmax.size() << ",";
                }
            }
            else
            {
                Sequence * seq = SequencesSetSPTree[j][i];
                if (seq->RC.empty() == false || (seq->Cmax.empty() == true && seq->RC.empty() == true)) {
                    fileRCsize << ",";
                    fileCmaxsize << ",";
                }
            }
        }
        fileRCsize << endl;
        fileCmaxsize << endl;
    }
}