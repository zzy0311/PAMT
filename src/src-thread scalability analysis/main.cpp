#include <iostream>
#include "ptmt.hpp"
#include "partition.hpp"
#include <omp.h>

int main(int argc, char * argv[]) {
    const auto t1 = chrono::steady_clock::now();
    
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <file> <Max_event> <Max_memory> <consecutive> <omega> <threads>" << std::endl;
        return 1;
    }

    string tmp (argv[1]);
    string gname = tmp.substr (tmp.find_last_of("/") + 1);
    int Max_event = stoi(argv[2]);
    int Max_memory = stoi(argv[3]); // delta
    string consecutive (argv[4]);
    int omega = stoi(argv[5]) * Max_memory;
    int threads = stoi(argv[6]);

    string out_file = "out_" + gname.substr(0, gname.size() - 4) + "_" + argv[2] + "_" + argv[3] + "_" + argv[4] + "_" + argv[5] + "_" + argv[6] + "_omp.txt";

    FILE* fp = fopen (out_file.c_str(), "w");

    // Set the number of threads for OpenMP
    omp_set_num_threads(threads);

    // Output configuration
    std::cout << "Max edge: " << Max_event << std::endl;
    std::cout << "delta: " << Max_memory << std::endl;
    std::cout << "Threads: " << threads << std::endl;

    // Read file and create a sorted list of temporal events
    vector<event> events; // event: {t, {u, v}}
    createEvents(tmp, events);

    std::cout << "# of events: " << events.size() << std::endl;
    const auto t2 = chrono::steady_clock::now();
    print_time(fp, "Read data time: ", t2 - t1);

    vector<tcons> ts_purn, ts_purn_repeat;
    int window_num = getTimeConstraint(events, ts_purn, ts_purn_repeat, (Max_event - 1) * Max_memory, omega);

    map<string, int> motif_count;

    // Parallel processing of main time windows
    #pragma omp parallel
    {
        map<string, int> local_motif_count;

        #pragma omp for schedule(dynamic)
        for (size_t j = 0; j < window_num; ++j) {
            vector<event> sub_events(events.begin() + ts_purn[j].first, 
                                     events.begin() + ts_purn[j].second + 1);

            instancemap sub_instances;
            set<vector<event>> local_keys;

            for (size_t i = 0; i < sub_events.size(); ++i) {
                countInstance(sub_events[i], sub_instances, local_keys, Max_event, Max_memory, consecutive);
            }

            for (auto it = sub_instances.begin(); it != sub_instances.end(); ++it) {
                string motif = encodeMotif(it->first);
                if (it->first.size() < 2) continue;
                local_motif_count[motif] += it->second.first;
            }
        }

        // Merge local results
        #pragma omp critical
        {
            for (const auto& pair : local_motif_count) {
                motif_count[pair.first] += pair.second;
            }
        }
    }

    // Parallel processing of repeated intervals
    #pragma omp parallel
    {
        map<string, int> local_motif_count;

        #pragma omp for schedule(dynamic)
        for (size_t j = 0; j < ts_purn_repeat.size(); ++j) {
            vector<event> sub_events(events.begin() + ts_purn_repeat[j].first, 
                                     events.begin() + ts_purn_repeat[j].second + 1);

            instancemap sub_instances;
            set<vector<event>> local_keys;

            for (size_t i = 0; i < sub_events.size(); ++i) {
                countInstance(sub_events[i], sub_instances, local_keys, Max_event, Max_memory, consecutive);
            }

            for (auto it = sub_instances.begin(); it != sub_instances.end(); ++it) {
                string motif = encodeMotif(it->first);
                if (it->first.size() < 2) continue;
                local_motif_count[motif] -= it->second.first;
            }
        }

        // Merge local results
        #pragma omp critical
        {
            for (const auto& pair : local_motif_count) {
                motif_count[pair.first] += pair.second;
            }
        }
    }

    const auto t3 = chrono::steady_clock::now();
    print_time(fp, "Count motifs time: ", t3 - t2);
    print_time(fp, "End-to-end Time: ", t3 - t1);

    for (auto it = motif_count.begin(); it != motif_count.end(); ++it) {
        fprintf(fp, "%s \t %d \n", it->first.c_str(), it->second);
    }

    fclose(fp);
    std::cout << "ALL DONE" << std::endl;
    return 0;
}
