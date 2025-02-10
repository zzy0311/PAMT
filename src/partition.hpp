#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "ptmt.hpp"

typedef pair<int, int> tcons;

int getTimeConstraint(vector<event>& events, 
                        vector<tcons>& ts_purn,   
                        vector<tcons>& ts_purn_repeat, 
                        int delta,
                        int omega);

// vector<event> getEventsSubset(vector<event>& events, int start, int end);

#endif
