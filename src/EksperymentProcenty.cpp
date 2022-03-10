

//
//int main() {
//
//    //string paths[] = {"..//crimes.txt"};
//    //string paths[] = {"..//crimes_modified_3.txt"};
//    //string paths[] = {"..//crimes_Pittsburgh.txt"};
//
//    //string path = paths[0];
//
//    //LoadDataset(path);
//
////
////    for (R = 350; R <= 350; R += 200) //I. 200, 300, 300 m,;; II. 500; 600 III. 100
////    {
////        for (T = 5760; T <=  5760; T += 14400) //1. 14400 (10 dni); 2. 11520 (8 dni); 3. 5760 (4 days);; II. 43200; 28800 III. 11520 (8 dni)
////        {
////            for (int i = 0; i < 1; i++) {
////                fstream exp3ExecTime;
////
////                double approxMargs[] = {0.001, 0.005, 0.01, 0.05, 0.1, 0.2};
////                int approxMargsSize = 6;
////
////                exp3ExecTime.open("..//Results//ExpPercentPittsburgh2 _R" + to_string(R) + "_T" + to_string(T) + "_"  + "crimes.txt",
////                                  fstream::out);
////
////                exp3ExecTime << "theta" << " ";
////
////                for (int k = 0; k < approxMargsSize; k++) {
////                    exp3ExecTime << approxMargs[k] << " ";
////                }
////
////                exp3ExecTime << "ST" << " " << "CS" << endl;
////
////                for (theta = 0.13; theta >= 0.045; theta -= 0.01) {
////
////                    exp3ExecTime << theta << " ";
////
////                    for (int k = 0; k < approxMargsSize; k++) {
////
////                        approxMargin = approxMargs[k];
////
////                        cout << "ACST_SPMiner_R" + to_string(R) + "_T" + to_string(T) + "_" + to_string(approxMargin) +
////                                "_" + paths[i] + ".txt" << endl;
////
////                        TransformData();
////                        SortDataset();
////
////                        ACST_SPMiner();
////                        int NumPatterns = CountPatternsACST_SPMiner();
////
////                        cout << theta << " " << NumPatterns << endl;
////                        exp3ExecTime << NumPatterns << " "; // for phase time
////
////                        ClearStructures();
////                    }
////
////
////                    cout << "CST_SPMiner_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;
////
////
////                    TransformData();
////                    SortDataset();
////
////
////                    CST_SPMiner();
////                    int NumPatterns = CountPatternsClosed();
////
////                    cout << theta << " " << NumPatterns << endl;
////
////
////                    exp3ExecTime << NumPatterns  << " ";
////                    ClearStructures();
////
////
////                    cout << "STBFM_R" + to_string(R) + "_T" + to_string(T) + "_" + paths[i] + ".txt" << endl;
////
////
////                    TransformData();
////                    SortDataset();
////
////
////                    STBFM();
////                    NumPatterns = CountPatterns();
////
////                    cout << " " << NumPatterns << endl;
////
////                    exp3ExecTime << NumPatterns  << " " << endl;
////                    ClearStructures();
////                }
////
////                exp3ExecTime.close();
////            }
////        }
////    }
//
//    for (R = 300; R <= 300; R += 200) //I. 200, 300, 300 m,;; II. 500; 600 III. 100
//    {
//        for (T = 11520; T <= 11520; T += 14400) //1. 14400 (10 dni); 2. 11520 (8 dni); 3. 5760 (4 days);; II. 43200; 28800 III. 11520 (8 dni)
//        {
//            for (int i = 0; i < 1; i++) {
//
//
//                double approxMarg = 0.05;
//
//                for (theta = 0.3; theta >= 0.3; theta -= 0.01) {
//
//
//                    cout << "ACST_SPMiner_R" + to_string(R) + "_T" + to_string(T) + "_" + to_string(approxMargin) +
//                            "_" + paths[i] + ".txt" << endl;
//
//                    TransformData();
//                    SortDataset();
//
//                    ACST_SPMiner();
//                    int NumPatterns = CountPatternsACST_SPMiner();
//                    PrintSequences();
//
//                    cout << theta << " " << NumPatterns << endl;
//
//
//                    ClearStructures();
//
//                }
//            }
//        }
//    }
//}
