//============================================================================
// Name        : aRticle.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <ctime>
#include <chrono>

#include "load_data.h"

using namespace std;
using namespace std::chrono;


//
// Created by Piotr on 03.11.2020.
//

///Experminets comp. no of patterns and time (Tables)
int main() {
    //1exp N 50,100 ... 250 crimes5000,crimes10000
    //2exp N 100, 200, 300, ... 1200, crimes10000
    //3exp theta 0.7 .. 0.2, 0.05 crimes3000


    //string paths[] = {"../crimes.txt"};
    //string paths[] = {"../crimes_modified_3.txt"};
    string paths[] = {"../crimes_Pittsburgh.txt"};

    string path = paths[0];

    LoadDataset(path);

    for (R = 350; R <= 350; R += 200) //I. 200, 300, 300 m,;; II. 500; 600 III. 100
    {
        for (T = 11520; T <= 11520; T += 14400) //1. 14400 (10 dni); 2. 11520 (8 dni); 3. 5760 (4 days);; II. 43200; 28800 III. 11520 (8 dni)
        {
            for (int i = 0; i < 1; i++) {
                fstream exp3ExecTime;

//
//                for(approxMargin = 0.025; approxMargin <= 0.11; approxMargin += 0.025) {
//
//                    exp3ExecTime.open("../results/CSTS_Miner_R" + to_string(R) + "_T"  + to_string(T) + "_" + to_string(approxMargin) + "_" + "_crimes.txt",
//                                      fstream::out);
//
//                    cout << "CSTS_Miner_R" + to_string(R) + "_T" + to_string(T) + "_" + to_string(approxMargin) + "_"  + paths[i] + ".txt" << endl;
//
//                    for (theta = 0.33; theta >= 0.245; theta -= 0.01) { //TimePhasesExp: 0.25, 0.255, 0.26, 0.265, 0.27
//                        //E1_pitssburgh theta = 0.33; theta >= 0.24; theta -= 0.01
//
//
//                        // string path = paths[i];
//                        // LoadDataset(path);
//                        TransformData();
//                        SortDataset();
//                        cout << "Loaded ACST-Miner" << endl;
//
//                        high_resolution_clock::time_point t1 = high_resolution_clock::now();
//                        CSTS_Miner();
//                        high_resolution_clock::time_point t2 = high_resolution_clock::now();
//
//                        int NumPatterns = CountPatternsCSTS_Miner();
//                        auto duration = duration_cast<seconds>(t2 - t1).count();
//                        cout << theta << " " << duration << " " << NumPatterns << endl;
//                        //cout << theta << " " << phase1Time << " " << phase2Time << " " << SequencesSetSPTree.size() << endl;
//
//                        //PrintSequences();
//
//                        exp3ExecTime << theta << " " << duration << " " << NumPatterns << endl; //for no of patterns
//                        //exp3ExecTime << theta << " " << phase1Time << " " << phase2Time << " " << SequencesSetSPTree.size() + 1 << endl; // for phase time
//                        ClearStructures();
//                    }
//
//                    exp3ExecTime.close();
//                }
//
//
//                //////////////
//
//                exp3ExecTime.open("../results/CTS_SPMiner_R" + to_string(R) + "_T" + to_string(T) + "_crimes.txt",
//                                  fstream::out);
//
//                cout << "CST_SPMiner_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;
//
//                for (theta = 0.33; theta >= 0.245; theta -= 0.01) {
//                    approxMargin = 0.4;
//
//                   // string path = paths[i];
//                   // LoadDataset(path);
//                    TransformData();
//                    SortDataset();
//                    cout << "Loaded CST-Miner" << endl;
//
//                    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//                    CST_SPMiner();
//                    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//
//                    int NumPatterns = CountPatternsClosed();
//                    auto duration = duration_cast<seconds>(t2 - t1).count();
//                    cout << theta << " " << duration << " " << NumPatterns << endl;
//
//                    PrintSequences();
//
//                    exp3ExecTime << theta << " " << duration << " " << NumPatterns << endl;
//                    ClearStructures();
//                }
//
//                exp3ExecTime.close();
//
//                //////////////
//
//                exp3ExecTime.open("../results/STBFM_R" + to_string(R) + "_T" + to_string(T) + "_crimes.txt",
//                                  fstream::out);
//
//                cout << "STBFM_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;
//
//                for (theta = 0.33; theta >= 0.245; theta -= 0.01) {
//                    approxMargin = 0.4;
//
//                   // string path = paths[i];
//                    //LoadDataset(path);
//                    TransformData();
//                    SortDataset();
//                    cout << "Loaded STBFM" << endl;
//
//                    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//                    STBFM();
//                    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//
//                    int NumPatterns = CountPatterns();
//                    auto duration = duration_cast<seconds>(t2 - t1).count();
//                    cout << theta << " " << duration << " " << NumPatterns << endl;
//
//                    PrintSequences();
//
//                    exp3ExecTime << theta << " " << duration << " " << NumPatterns << endl;
//                    ClearStructures();
//                }
//
//                exp3ExecTime.close();


//                //////////////
//
//                exp3ExecTime.open("../results/CSTPM_R" + to_string(R) + "_T" + to_string(T) + "_crimes.txt",
//                                  fstream::out);
//
//                cout << "CSTPM_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;
//
//               for (theta = 0.33; theta >= 0.245; theta -= 0.01) {
//                    approxMargin = 0.4;
//
//                    // string path = paths[i];
//                    //LoadDataset(path);
//                    TransformData();
//                    SortDataset();
//                    cout << "Loaded CSTPM" << endl;
//
//                    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//                    CSTPM();
//                    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//
//                    int NumPatterns = CountPatterns();
//                    auto duration = duration_cast<seconds>(t2 - t1).count();
//                    cout << theta << " " << duration << " " << NumPatterns << endl;
//
//                    PrintSequences();
//
//                    exp3ExecTime << theta << " " << duration << " " << NumPatterns << endl;
//                    ClearStructures();
//                }
//
//                exp3ExecTime.close();

                //////////////

                exp3ExecTime.open("../results/STS-Miner_R" + to_string(R) + "_T" + to_string(T) + "_crimes.txt",
                                  fstream::out);

                cout << "STS-Miner_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;

                for (theta = 0.33; theta >= 0.245; theta -= 0.01) {
                    approxMargin = 0.4;

                    // string path = paths[i];
                    //LoadDataset(path);
                    TransformData();
                    SortDataset();
                    cout << "Loaded STS-Miner" << endl;

                    high_resolution_clock::time_point t1 = high_resolution_clock::now();
                    STS_Miner();
                    high_resolution_clock::time_point t2 = high_resolution_clock::now();

                    int NumPatterns = CountPatterns();
                    auto duration = duration_cast<seconds>(t2 - t1).count();
                    cout << theta << " " << duration << " " << NumPatterns << endl;

                    PrintSequences();

                    exp3ExecTime << theta << " " << duration << " " << NumPatterns << endl;
                    ClearStructures();
                }

                exp3ExecTime.close();

            }
        }
    }


}