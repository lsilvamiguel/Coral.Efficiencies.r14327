#ifndef CsPPI_EC02tools_H
#define CsPPI_EC02tools_H

int ReadFile(const char* fname,const bool first,int& global_run_nmb,unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]);

void init_d_array(unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]);
int WriteFile(const char* fname,int &global_run_nmb,unsigned int data_array[200][6][8][8][5],unsigned int global_checksum[200]);
#endif
