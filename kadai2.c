#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 //与えられるプロモータ領域の最大遺伝子数
#define ACGT 4
#define threshold 6 //閾値

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

int frequency_matrix(int transcription_len, float frequency[ACGT][transcription_len]){
  int i=0, j=0;
  for(i=0; i<ACGT; i++){
    for(j=0; j<transcription_len; j++){
        frequency[i][j]=0;
    }
  }

  for(j=0; j<transcription_len; j++){
    for(i=0; i<MAX_SEQ_NUM; i++){
      switch(g_motif[i][j]){
        case 'A' : frequency[0][j]++; break;
        case 'C' : frequency[1][j]++; break;
        case 'G' : frequency[2][j]++; break;
        case 'T' : frequency[3][j]++; break;
      }
    }
  }
}

int score_matrix(int transcription_len, float score[ACGT][transcription_len]){
  int i=0, j=0;
  float m = score[0][0]+score[1][0]+score[2][0]+score[3][0]+4;
  for(i=0; i<ACGT; i++){
    for(j=0; j<transcription_len; j++){
        score[i][j]=score[i][j]+1;
    }
  }
  for(i=0; i<ACGT; i++){
    for(j=0; j<transcription_len; j++){
        score[i][j]=score[i][j]/m;
    }
  }

  float qAT = 7519429.0/(7519429+4637676+4637676+7519429);
  float qCG = 4637676.0/(7519429+4637676+4637676+7519429);
  float q[ACGT]={qAT, qCG, qCG, qAT};

  for(i=0; i<ACGT; i++){
    for(j=0; j<transcription_len; j++){
        score[i][j]=log(score[i][j]/q[i]);
    }
  }
}

int binding_site(int transcription_len,  float binding[ACGT][transcription_len]){
  int i=0, j=0, k=0;
  int promoter_len=strlen(g_pro[0].seq);
  int pos=0;
  float score=0.0;
  char array[transcription_len+1];
  for(i=0; i<MAX_GENE_NUM; i++){
    for(j=0; j<promoter_len-transcription_len; j++){
      score=0.0;
      for(k=0; k<transcription_len; k++){
        array[k]=g_pro[i].seq[j+k];
        switch(g_pro[i].seq[j+k]){
          case 'A' : score+=binding[0][k]; break;
          case 'C' : score+=binding[1][k]; break;
          case 'G' : score+=binding[2][k]; break;
          case 'T' : score+=binding[3][k]; break;
        }
      }
      array[transcription_len] = '\0';
      
      if (score >= threshold) {
        printf("pro: %s\n",g_pro[i].name);
        printf("pos: %d\n",j+1);
        printf("hit(%s)= %f\n\n",array,score);
      }
    }
  }
  
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む

  int transcription_factor_len=strlen(g_motif[0]);
  float matrix[ACGT][transcription_factor_len];

  frequency_matrix(transcription_factor_len, matrix);
  score_matrix(transcription_factor_len, matrix);
  binding_site(transcription_factor_len, matrix);
  
  return 0;
}