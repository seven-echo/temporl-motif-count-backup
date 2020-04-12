#ifndef snap_temporalmotifsampler_h
#define snap_temporalmotifsampler_h

#include "Snap.h"

#include "temporalmotiftypes.h"

#include <vector>
#include <map>
#include <stack>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <limits>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>

// TODO: no global namespaces
using namespace std;

// TODO: Move this to cpp file
class MotifSampler {
public:
	MotifSampler() {}

	TUInt64 CountMotifs() {
		return 0;
	}
};

// TODO: REPLACE THIS WITH SNAP STRUCTURES
struct TEdge {
	long long src, dst, tim, id;
	const bool operator<(const TEdge& o) const {
		if (o.tim != tim) return tim < o.tim;
		if (o.src != src) return src < o.src;
		if (o.dst != dst) return dst < o.dst;
		if (o.id != id) return id < o.id;
		return false;
	}

	friend ostream& operator<<(ostream& in, const TEdge& o) {
		in << "(" << o.src << "," << o.dst << "," << o.tim << "," << o.id << ")";
		return in;
	}
};

class TempMotifSampler {
public:
    MotifSampler ms;
    TempMotifSampler(const TStr& filenameG, const TStr& filenameM);

//    double ExactEnumerateMotifs(int delta, int window, vector<TEdge>& edges) {
    std::tuple<double, map<pair<int, int>, double>> ExactEnumerateMotifs(int delta, int window, vector<TEdge>& edges) {
    	int nV = adj_list_.size();
    	edgeCount_.clear();
    	edgeCount_.resize(nV, 0);
    	mapGM_.clear();
    	mapGM_.resize(nV, -1);
    	mapMG_.clear();
    	mapMG_.resize(Vm_, -1);

    	double res = 0;
    	vector<int> eStack;

    	int eG = 0, eM = 0, uG = -1, vG = -1, uM = -1, vM = -1;
    	const int INF = numeric_limits<int>::max();
    	int _t = INF;
	// +++
	map<pair<int, int>, double> node_counter;
	// +++
    	while (true) {
    		int last = eStack.empty() ? -1 : eStack.back();
    		eG = FindNextMatch(eM, eG, mapMG_, mapGM_, _t, edges, last);
    		TEdge edge = {0, 0, 0, INF};
    		if (eG < (int) edges.size()) {
    			if (eM == int(edgesM_.size()) - 1) {
			//if (eM == int(edgesM_.size())) {
    				// Apply the weight function
    				if (window == -1) {
					res += 1;
					for (auto e : eStack){
						TEdge temp_edge = edges[e];
						uG = temp_edge.src, vG = temp_edge.dst;
						//cout << uG << "\t" << vG << endl;
						node_counter[make_pair(uG, vG)] += 1;
						//cout << e << endl;
					}
					// Add last edge
					TEdge temp_edge = edges[eG];
					uG = temp_edge.src, vG = temp_edge.dst;
					//cout << uG << "\t" << vG << endl;
					node_counter[make_pair(uG, vG)] += 1;
					//cout << eG <<endl;
					//cout << "match once!" << endl;
				}
    				else {
					res += 1. / (1 - double(edges[eStack.back()].tim - edges[eStack[0]].tim) / window);
					for (auto e : eStack){
						TEdge temp_edge = edges[e];
						uG = temp_edge.src, vG = temp_edge.dst;
						node_counter[make_pair(uG, vG)] += 1. / (1 - double(edges[eStack.back()].tim - edges[eStack[0]].tim) / window);
						// cout << e << endl;
					}
				};
    			} else {
    				edge = edges[eG];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM_[eM].first, vM = edgesM_[eM].second;

    				mapGM_[uG] = uM;
    				mapGM_[vG] = vM;
    				mapMG_[uM] = uG;
    				mapMG_[vM] = vG;

    				edgeCount_[uG] += 1;
    				edgeCount_[vG] += 1;
    				if (eStack.empty()) {
    					_t = edge.tim + delta;
    				}
    				eStack.push_back(eG);
    				eM += 1;
				//cout << eM << endl;
    			}
    		}
    		eG += 1;

    		while (eG >= (int) edges.size() || edge.tim > _t) {
    			if (!eStack.empty()) {
    				eG = eStack.back() + 1;
    				eStack.pop_back();

    				edge = edges[eG - 1];
    				uG = edge.src, vG = edge.dst;
    				uM = edgesM_[eM].first, vM = edgesM_[eM].second;

    				if (eStack.empty()) {
    					_t = INF;
    				}
    				edgeCount_[uG] -= 1;
    				edgeCount_[vG] -= 1;

    				if (edgeCount_[uG] == 0) {
    					uM = mapGM_[uG];
    					mapMG_[uM] = -1;
    					mapGM_[uG] = -1;
    				}

    				if (edgeCount_[vG] == 0) {
    					vM = mapGM_[vG];
    					mapMG_[vM] = -1;
    					mapGM_[vG] = -1;
    				}

    				eM -= 1;
    			} else {

    				//return res;
// +++
				return std::make_tuple(res, node_counter);
// +++
    			}
    		}
    	}
    }

    inline int FindNextMatch(
    	int eM, 
    	int eG, 
    	vector<int>& mapMG, 
    	vector<int>& mapGM, 
    	int _t,
    	vector<TEdge>& edges,
    	int last) {

    	int uM, vM, uG, vG;
    	uM = edgesM_[eM].first;
    	vM = edgesM_[eM].second;
    	uG = mapMG_[uM];
    	vG = mapMG_[vM];

    	vector<TEdge>* S;
    	int head = 0;
    	if (uG >= 0 && vG >= 0) {
    		S = &adjMap_[uG][vG];
    	} else if (uG >= 0) {
    		S = &adj_list_[uG];
    	} else if (vG >= 0) {
    		S = &revadj_list_[vG];
    	} else {
    		S = &edges;
    	}
    	bool small = (S->size() < 16);
    	if (!small) {
	    	head = lower_bound(S->begin(), S->end(), edges[eG]) - S->begin();
	    }

    	for (int i = head; i < (int) S->size(); i++) {
    		auto edge = (*S)[i];
    		if (edge.id < eG || edge.tim > _t) {
    			if (small) continue;
    			else break;
    		}
    		if (last != -1 && edge.tim <= edges[last].tim) {
    			continue;
    		}

    		int _eG = edge.id;
    		int _uG = edge.src, _vG = edge.dst;
    		if (uG == _uG || (uG < 0 && mapGM_[_uG] < 0)) {
    			if (vG == _vG || (vG < 0 && mapGM_[_vG] < 0)) {
    				return _eG;
    			}
    		} 
    	}
    	return edges.size();
    }

//    TUInt64 ExactCountMotifs(int delta) {
    std::tuple<double, map<pair<int, int>, double>> ExactCountMotifs(int delta) {
    	InitializeStructures(edges_);
	
    	return ExactEnumerateMotifs(delta, -1, edges_);
    }

    double ApproximateCountMotifsRandomWindow(int delta) {
    	srand(time(nullptr));

    	// Randomly partition graph
		const int window = 100 * delta;
		int n_trials = 200;
		
		double tot_estimate = 0;
		map<pair<int, int>, double> counter;
		vector<TEdge> edges;
		for (int c = 0; c < n_trials; c++) {
			edges.clear();

			int offset = (rand() % window) - window / 2;

			// pick a random edge and a random delta window around it
			int ei = rand() % edges_.size();

			int ti = edges_[ei].tim;
			int st = ((ti-offset > 0 ? ti-offset : ti - window + 1 - offset)/window) * window + offset + window;

			TEdge curr = {0, 0, st - window, 0};
			int idx = lower_bound(edges_.begin(), edges_.end(), curr) - edges_.begin();
			while (idx < (int) edges_.size() && edges_[idx].tim <= st) {
				edges.push_back(edges_[idx]);
				idx++;
			}
			InitializeStructures(edges);

			double pi = edges.size() / double(edges_.size());
// +++
			auto my_temp = ExactEnumerateMotifs(delta, window, edges);
			double res = std::get<0>(my_temp) / pi;
			auto node_counter = std::get<1>(my_temp);
			for(auto& t_counter : node_counter){
          			//node_counter[t_counter.first] = t_counter.second / pi;
				// update counter
				counter[t_counter.first] += t_counter.second / pi;
        		} 
// +++
//			double res = ExactEnumerateMotifs(delta, window, edges) / pi;

			tot_estimate += res;
		}
//TODO +++ save counter file
	// +++
	printf("save file......\n");
	FILE *fp=fopen("./motif_1.txt","w");
	// 检测文件是否打开成功；
	if(!fp){
		printf("open file failed\n");
		return -1; //返回异常
	}
	for(auto& t_counter : counter){
		// cout << t_counter.first.first;
		// printf("%d",t_counter.first.first);
		fprintf(fp, "%d %d %.f\n", t_counter.first.first, t_counter.first.second, t_counter.second / n_trials);
	}
	fclose(fp);
	//+++
//
    	return tot_estimate / n_trials;
    }

    double ApproximateCountMotifsSlidingWindowSkip(int delta, double c = 3e1, double mult = 3e1) {
		srand(time(nullptr));

		default_random_engine generator;
		uniform_real_distribution<double> distribution(0.0, 1.0);
    	
		const int window = c * delta;
		const int n_trials = 1;

		double tot_estimate = 0;
		map<pair<int, int>, double> counter;
		vector<TEdge> edges;
		for (int c = 0; c < n_trials; c++) {
			edges.clear();
			int offset = rand() % window;

			double res = 0;
			map<pair<int, int>, double> node_counter;
			int nxt = edges_[0].tim + offset, idxL = 0, idxR = 0;
			while (idxR < (int) edges_.size() && edges_[idxR++].tim <= nxt);
			while (idxL < (int) edges_.size()) {
				int nedges = idxR - idxL;

				double pi = min(mult * double(nedges) / edges_.size(), 1.0);
				double p = distribution(generator);
				if (p <= pi) {
					edges.assign(edges_.begin() + idxL, edges_.begin() + idxR);
					InitializeStructures(edges);
// +++
					auto my_temp = ExactEnumerateMotifs(delta, window, edges);
					auto res = std::get<0>(my_temp) / pi;
					node_counter = std::get<1>(my_temp);
					for(auto& t_counter : node_counter){
			  			node_counter[t_counter.first] = t_counter.second / pi;
					} 
// +++
					//res += ExactEnumerateMotifs(delta, window, edges) / pi;
					edges.clear();
				}

				idxL = idxR;
				if (idxR < (int) edges_.size()) {
					int t = edges_[idxR].tim;
					nxt = ((t-offset)/window) * window + window + offset;
				}
				while (idxR < (int) edges_.size() && edges_[idxR++].tim <= nxt);
			}
			tot_estimate += res;
//+++
			for(auto& t_counter : node_counter){
			  	counter[t_counter.first] += t_counter.second;
			}
//+++
		}
//TODO save counter to file
	// +++
	printf("save file......\n");
	FILE *fp=fopen("motif_1_.txt","w");
	// 检测文件是否打开成功；
	if(!fp){
		printf("open file failed\n");
		return -1; //返回异常
	}
	for(auto& t_counter : counter){
		// cout << t_counter.first.first;
		// printf("%d",t_counter.first.first);
		fprintf(fp, "%d %d %.f\n", t_counter.first.first, t_counter.first.second, t_counter.second / n_trials);
	}
	fclose(fp);
	//+++
//
    	return tot_estimate / n_trials;
    }

    double ApproximateCountMotifs(int delta) {
    	return ApproximateCountMotifsSlidingWindowSkip(delta);
    }
    
private:
	void InitializeStructures(vector< TEdge >& edgeList) {
		unordered_map<int, int> remap;
		int id = 0;
		for (auto e : edgeList) {
			if (!remap.count(e.src)) {
				remap[e.src] = id++;
			}
			if (!remap.count(e.dst)) {
				remap[e.dst] = id++;
			}
		}
		for (auto& e : edgeList) {
			e.src = remap[e.src];
			e.dst = remap[e.dst];
		}
		sort(edgeList.begin(), edgeList.end());

		for (int i = 0; i < (int) edgeList.size(); i++) {
			edgeList[i].id = i;
		}

		// TODO: IMPLEMENT SAMPLING MORE ELEGANTLY
		adj_list_.clear();
		revadj_list_.clear();
		adjMap_.clear();

		adj_list_.resize(id);
		revadj_list_.resize(id);
		adjMap_.resize(id);
		for (const auto& edge : edgeList) {
			adj_list_[edge.src].push_back(edge);
			revadj_list_[edge.dst].push_back(edge);
			adjMap_[edge.src][edge.dst].push_back(edge);
		}
	} 

	// Edge index based adjacency list
	vector< vector< TEdge > > adj_list_;
	vector< vector< TEdge > > revadj_list_;
	vector< unordered_map<int, vector<TEdge> > > adjMap_;

	// Structures for algo
	vector<int> edgeCount_;
	vector<int> mapGM_;
	vector<int> mapMG_;

	// TODO: REPLACE WITH SNAP STRUCTURES
	vector< TEdge > edges_;
	int Vm_;
	vector< pair<int, int> > edgesM_;
};

#endif  // snap_temporalmotifsampler_h
