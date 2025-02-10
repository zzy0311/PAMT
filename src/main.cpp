#include <iostream>
#include "ptmt.hpp"
#include "partition.hpp"
#include <omp.h>

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    int Max_event = stoi(argv[2]);
    int Max_memory = stoi(argv[3]); // delta
    
    string consecutive (argv[4]);
    string out_file = "out_" + gname.substr(0,gname.size()-4) + "_" + argv[2] + "_" + argv[3] + "_" + argv[4] + "_"+argv[5]+"_omp.txt";

    FILE* fp = fopen (out_file.c_str(), "w");
    
    cout << "Max edge: " << Max_event << endl;
    
    cout << "delta: " << Max_memory << endl;
    
//  Read file and create a sorted list of temporal events
    vector<event> events; //event: {t, {u, v}}
    createEvents(tmp, events);
        
    cout << "# of events:" << events.size() << endl;
    const auto t2 = chrono::steady_clock::now();
    print_time (fp, "Read data time: ", t2 - t1);



//  #####################################################
    int omega = stoi(argv[5]) ;//omega
    vector<tcons> ts_purn, ts_purn_repeat;                                 // l*delta              //omega
    int window_num = getTimeConstraint(events, ts_purn, ts_purn_repeat, (Max_event) * Max_memory, omega);
 const auto t3 = chrono::steady_clock::now();
    map<string, int> motif_count;

   
    #pragma omp parallel
    {
        map<string, int> local_motif_count;
        
        #pragma omp for schedule(dynamic)
        for(size_t j = 0; j < window_num; ++j) {
            vector<event> sub_events(events.begin() + ts_purn[j].first, 
                                   events.begin() + ts_purn[j].second + 1);
            
            instancemap sub_instances;
            set<vector<event>> local_keys;
            
            for(size_t i = 0; i < sub_events.size(); ++i) {
                countInstance(sub_events[i], sub_instances, local_keys,
                            Max_event, Max_memory, consecutive);
            }

            for(auto it = sub_instances.begin(); it != sub_instances.end(); ++it) {
                string motif = encodeMotif(it->first);
                if (it->first.size() < 2) continue;
                local_motif_count[motif] += it->second.first;
            }
        }

        
        #pragma omp critical
        {
            for(const auto& pair : local_motif_count) {
                motif_count[pair.first] += pair.second;
            }
        }
    }

    
    #pragma omp parallel
    {
        map<string, int> local_motif_count;
        
        #pragma omp for schedule(dynamic)
        for(size_t j = 0; j < ts_purn_repeat.size(); ++j) {
            vector<event> sub_events(events.begin() + ts_purn_repeat[j].first, 
                                   events.begin() + ts_purn_repeat[j].second + 1);
            
            instancemap sub_instances;
            set<vector<event>> local_keys;
            
            for(size_t i = 0; i < sub_events.size(); ++i) {
                countInstance(sub_events[i], sub_instances, local_keys,
                            Max_event, Max_memory, consecutive);
            }

            for(auto it = sub_instances.begin(); it != sub_instances.end(); ++it) {
                string motif = encodeMotif(it->first);
                if (it->first.size() < 2) continue;
                local_motif_count[motif] -= it->second.first;
            }
        }

      
        #pragma omp critical
        {
            for(const auto& pair : local_motif_count) {
                motif_count[pair.first] += pair.second;
            }
        }
    }

//  ######################################################



    
    const auto t4 = chrono::steady_clock::now();
    print_time (fp, "Count motifs time: ", t4 - t3);
    print_time (fp, "partition Time: ", t3 - t2);
    
    for (auto it=motif_count.begin(); it!=motif_count.end(); ++it) {
        fprintf(fp, "%s \t %d \n", (*it).first.c_str(), (*it).second);
    }
    
    fclose (fp);
    cout << "ALL DONE" << endl;
    return 0;
}
