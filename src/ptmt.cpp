#include "ptmt.hpp"
#include <omp.h>

void countInstanceParallel(event e, instancemap& imap, set<vector<event>>& keys, 
                          int N_event, int d_c, string consecutive,
                          omp_lock_t& map_lock) {
    bool ie {true};
    vertex u = e.second.first;
    vertex v = e.second.second;
    vector<vector<event>> new_motif;

    
    for (auto it = keys.begin(); it != keys.end();) {
        vector<event> key = *it;
        if (e.first - key.back().first <= d_c) {
            if (key.size() < N_event) {
                set<vertex> nodes;
                
                
                omp_set_lock(&map_lock);
                nodes = imap[key].second;
                omp_unset_lock(&map_lock);
                
                nodes.insert(u);
                nodes.insert(v);
                if (nodes.find(u) != nodes.end() || nodes.find(v) != nodes.end()) {
                    if (key.back().first != e.first || true) {
                        vector<event> motif = key;
                        motif.push_back(e);
                        new_motif.push_back(motif);
                        
                        omp_set_lock(&map_lock);
                        imap[motif].first += imap[key].first;
                        imap[motif].second = nodes;
                        omp_unset_lock(&map_lock);
                        
                        ie = false;
                        if (consecutive == "YES") {
                            it = keys.erase(it);
                            continue;
                        }
                    }
                }
                ++it;
            } else {
                it = keys.erase(it);
            }
        } else {
            it = keys.erase(it);
        }
    }

    if (!new_motif.empty()) {
        omp_set_lock(&map_lock);
        for (vector<event> const &mt: new_motif) {
            keys.insert(mt);
        }
        omp_unset_lock(&map_lock);
    }

    vector<event> E;
    E.push_back(e);
    
    omp_set_lock(&map_lock);
    imap[E].first += 1;
    imap[E].second.insert(u);
    imap[E].second.insert(v);
    keys.insert(E);
    omp_unset_lock(&map_lock);
}

void createEvents (string filename, vector<event>& events){
    ifstream in(filename);
    string line;
    
    while (getline(in, line)) {
        if (line[0] != '%' && line[0] != '#'){
            stringstream ss (line);
            vertex u, v;
            timestamp t;
            edge e;
            ss >> u >> v >> t;
            if (u != v) {
                e = make_pair(u, v);
                events.push_back(make_pair(t, e));
            }
        }
    }
    sort(events.begin(), events.end());
    events.erase(unique(events.begin(), events.end()),events.end()); //remove duplicates
    return;
}

void countInstance (event e, instancemap& imap, set<vector<event>>& keys, int N_event, int d_c, string consecutive){
    bool ie {true};
    vertex u = e.second.first;
    vertex v = e.second.second;
     //#typedef unordered_map<vector<event>, pair<int, set<vertex>>> instancemap; //a hashtable of key and counts
    vector<vector<event>> new_motif;    //used to store the new motifs
      //#typedef  set<vector<event>> keys;
    for (auto it = keys.begin(); it != keys.end();) {   //for each current prefix
        vector<event> key = *it;
        if (e.first - key.back().first <= d_c) {  //check delta C and delta W
            if (key.size() < N_event) { //check the number of events
                set<vertex> nodes = imap[key].second;
                nodes.insert(u);
                nodes.insert(v);
                if (imap[key].second.find(u)!=imap[key].second.end() || imap[key].second.find(v)!=imap[key].second.end()) {
                    if (key.back().first!=e.first||true) { //check synchronous events
                        vector<event> motif = key;
                        motif.push_back(e);
                        new_motif.push_back(motif);
                        imap[motif].first += imap[key].first;
                        imap[motif].second = nodes;
                        ie = false;
                        if (consecutive == "YES") {
                            it = keys.erase(it);
                            continue;
                        }
                    }
                }
                ++it;
            } else {
                it = keys.erase(it);    //remove prefix if it exceeds the size constrain
            }
        } else {
            it = keys.erase(it);    //remove prefix if it exceeds the delta constrain
        }
    }
    //add the new motifs to the current prefix list
    if (!new_motif.empty()) {
        for (vector<event> const &mt: new_motif) {
            keys.insert(mt);
        }
    }

    vector<event> E;
    E.push_back(e);
    imap[E].first += 1;
    imap[E].second.insert(u);
    imap[E].second.insert(v);
    keys.insert(E); // add the new event to the current prefix list
    
    return;
}


string encodeMotif(vector<event> instance){
    string motif;
    string temp;
    map<vertex, string> code;
    int i=0;
    for (auto it=instance.begin(); it!=instance.end(); ++it) {
        vertex u = it->second.first;
        vertex v = it->second.second;
        timestamp t = it->first;
        if (code.find(u)==code.end()) {
            code[u] = to_string(i);
            i++;
        }
        temp.append(code[u]);
        if (code.find(v)==code.end()) {
            code[v] = to_string(i);
            i++;
        }
        temp.append(code[v]);
    }
    motif.append(temp);

    return motif;
}


char sconvert (int i) {
    string s("abcdefghijklmnopqrstuvwxyz");
    return s.at(i-1);
}

set<vertex> getNodes(vector<event> key){
    set<vertex> nodes;
    for (int i=0; i<key.size(); i++) {
        nodes.insert(key[i].second.first);
        nodes.insert(key[i].second.second);
    }
    return nodes;
}

