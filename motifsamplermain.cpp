#include "temporalmotifs.h"
#include "temporalmotifsampler.h"

#ifdef USE_OPENMP
#include <omp.h>
#include <string.h>
#endif

#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

int main(int argc, char* argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);

  const TStr temporal_graph_filename =
    Env.GetIfArgPrefixStr("-i:", "simple-example.txt",
			  "Input directed temporal graph file");
  const TStr motif_graph_filename =
    Env.GetIfArgPrefixStr("-m:", "simple-motif.txt",
        "Input directed motif graph file");
  const TStr output =
    Env.GetIfArgPrefixStr("-o:", "motif-count.txt",
        "Output motif count file");
  const TFlt delta =
    Env.GetIfArgPrefixFlt("-delta:", 4096, "Time window delta");
  const TFlt c =
    Env.GetIfArgPrefixFlt("-c:", 30, "Window size multiplier");
  const TFlt mult =
    Env.GetIfArgPrefixFlt("-r:", 30, "Sampling probability multiplier");

  TempMotifSampler tmc(temporal_graph_filename, motif_graph_filename);

  Env.PrepArgs(TStr::Fmt("Temporalmotifs. build: %s, %s. Time: %s",
       __TIME__, __DATE__, TExeTm::GetCurTm()));  

  double curr = get_wall_time(); 
//
  printf("ExactCountMotifs sum motifs: %d\n", (int) std::get<0>(tmc.ExactCountMotifs(delta)));
  // +++
  auto counter = std::get<1>(tmc.ExactCountMotifs(delta));
  printf("save file......\n");

  FILE *fp=fopen(output.CStr(),"w");
  // 检测文件是否打开成功；
  if(!fp){
	printf("open file failed\n");
	return -1; //返回异常
  }
  for(auto& t_counter : counter){
	// cout << t_counter.first.first;
	// printf("%d",t_counter.first.first);
	fprintf(fp, "%d %d %.f\n", t_counter.first.first, t_counter.first.second, t_counter.second);
  }
  fclose(fp);
//+++
//
//
  //printf("%d\n", (int) tmc.ExactCountMotifs(delta));
  //printf("%d\n", (int) tmc.BruteCountM2(delta));
  //printf("%d\n", (int) tmc.BruteCountM2ReallySlow(delta));
  printf("\nrun time: %lfs\n", get_wall_time() - curr);

/*
  TExeTm ExeTm2;
  printf("%lf\n", tmc.ApproximateCountMotifsSlidingWindow(delta));
  printf("\nrun time: %s (%s)\n", ExeTm2.GetTmStr(),
   TSecTm::GetCurTm().GetTmStr().CStr());
//*/

/*
  TExeTm ExeTm3;
  printf("%lf\n", tmc.ApproximateCountMotifsRandomWindow(delta));
  printf("\nrun time: %s (%s)\n", ExeTm3.GetTmStr(),
	 TSecTm::GetCurTm().GetTmStr().CStr());
//*/

//  curr = get_wall_time();
//  printf("ApproximateCountMotifsSlidingWindowSkip motif sum number:\t%lf\n", tmc.ApproximateCountMotifsSlidingWindowSkip(delta, c, mult));
//  printf("\nrun time: %lfs\n", get_wall_time() - curr);
  return 0;
}
