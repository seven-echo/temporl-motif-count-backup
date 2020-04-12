#include "Snap.h"

#include "temporalmotifsampler.h"

#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
TempMotifSampler::TempMotifSampler(const TStr& filenameG, const TStr& filenameM) {
  // Formulate input File Format:
  //   source_node destination_node timestamp
  TTableContext context;
  Schema temp_graph_schema;
  temp_graph_schema.Add(TPair<TStr,TAttrType>("source", atInt));
  temp_graph_schema.Add(TPair<TStr,TAttrType>("destination", atInt));
  temp_graph_schema.Add(TPair<TStr,TAttrType>("time", atInt));

  // Load the temporal graph
  PTable data_ptr = TTable::LoadSS(temp_graph_schema, filenameG, &context, ' ');
  TInt src_idx = data_ptr->GetColIdx("source");
  TInt dst_idx = data_ptr->GetColIdx("destination");
  TInt tim_idx = data_ptr->GetColIdx("time");
  set<long long> input_times;
  set<long long> nodes;
  for (TRowIterator RI = data_ptr->BegRI(); RI < data_ptr->EndRI(); RI++) {
    TInt row_idx = RI.GetRowIdx();
    long long src = data_ptr->GetIntValAtRowIdx(src_idx, row_idx).Val;
    long long dst = data_ptr->GetIntValAtRowIdx(dst_idx, row_idx).Val;
    long long tim = data_ptr->GetIntValAtRowIdx(tim_idx, row_idx).Val;
    // Do not include self loops as they do not appear in the definition of
    // temporal motifs.
    if (src != dst) { 
      edges_.push_back({src, dst, tim, 0LL});
      input_times.insert(tim);
      nodes.insert(src);
      nodes.insert(dst);
    }
  }
  cerr << "Max time is: " << *(--input_times.end()) << endl;
  cerr << "Max node is: " << *(--nodes.end()) << endl;
  if (input_times.size() != edges_.size()) {
    cerr << "Discreptancy size: " << input_times.size() << " " << edges_.size() << endl;
    //throw exception();//"Error! Input time stamps are not unique!");
  }
  sort(edges_.begin(), edges_.end());

  for (int i = 0; i < (int) edges_.size(); i++) {
    //cerr << edges_[i] << " became " << i << endl;
    edges_[i].id = i;
  }

  InitializeStructures(edges_);

  // Load the motif graph
  int max_nodes = 0;
  data_ptr = TTable::LoadSS(temp_graph_schema, filenameM, &context, ' ');
  src_idx = data_ptr->GetColIdx("source");
  dst_idx = data_ptr->GetColIdx("destination");
  tim_idx = data_ptr->GetColIdx("time");
  for (TRowIterator RI = data_ptr->BegRI(); RI < data_ptr->EndRI(); RI++) {
    TInt row_idx = RI.GetRowIdx();
    int src = data_ptr->GetIntValAtRowIdx(src_idx, row_idx).Val;
    int dst = data_ptr->GetIntValAtRowIdx(dst_idx, row_idx).Val;
    // Do not include self loops as they do not appear in the definition of
    // temporal motifs.
    if (src != dst) { 
      edgesM_.push_back(make_pair(src, dst));
    }
    max_nodes = max(max_nodes, src + 1);
    max_nodes = max(max_nodes, dst + 1);
  }
  Vm_ = max_nodes;
}
