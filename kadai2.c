#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define NUM 4

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }
  
  int num = strlen(g_motif[0]);
  float hindo[4][num];
  int i=0, j=0;
  for(i=0; i<4; i++){
    for(j=0; j<num; j++){
        hindo[i][j]=0;
    }
  }
  for(j=0; j<num; j++){
    for(i=0; i<MAX_SEQ_NUM; i++){
        switch(g_motif[i][j]){
            case 'A' : hindo[0][j]++; break;
            case 'C' : hindo[1][j]++; break;
            case 'G' : hindo[2][j]++; break;
            case 'T' : hindo[3][j]++; break;
        }
    }
  }

  /*for(i=0; i<4; i++){
  for(j=0; j<num; j++){
    printf("%d ",hindo[i][j]+1);
    }
    printf("\n");
  }*/

  float m = hindo[0][0]+hindo[1][0]+hindo[2][0]+hindo[3][0]+4;
  for(i=0; i<4; i++){
    for(j=0; j<num; j++){
        hindo[i][j]=hindo[i][j]+1;
    }
  }
  for(i=0; i<4; i++){
    for(j=0; j<num; j++){
        hindo[i][j]=hindo[i][j]/m;
    }
  }
  

/*for(i=0; i<4; i++){
  for(j=0; j<num; j++){
    printf("%f ",hindo[i][j]);
    }
    printf("\n");
  }*/

  float qAT = 7519429.0/(7519429+4637676+4637676+7519429);
  float qCG = 4637676.0/(7519429+4637676+4637676+7519429);
  float q[4]={qAT, qCG, qCG, qAT};

  float hindo_1[4][num];
  for(i=0; i<4; i++){
    for(j=0; j<num; j++){
        hindo_1[i][j]=0;
    }
  }
  for(i=0; i<4; i++){
    for(j=0; j<num; j++){
        hindo_1[i][j]=log(hindo[i][j]/q[i]);
    }
  }

  for(i=0; i<4; i++){
  for(j=0; j<num; j++){
    printf("%f ",hindo_1[i][j]);
    }
    printf("\n");
  }

  int len = strlen(g_pro[0].seq);
  

  return 0;
}
