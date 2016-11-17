///////////////////////////////////////////////////////////////////////////////
// read txt file , make ntuple out of it
///////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include "TObject.h"
#include "TNtuple.h"

//-----------------------------------------------------------------------------
void parse_number(const char* Token, float& Var) {
  if ((index(Token,'e') < 0) && (index(Token,'E') < 0)) sscanf(Token,"%f",&Var);
  else                                                  sscanf(Token,"%e",&Var);
}

//-----------------------------------------------------------------------------
int create_ntuple(const char* Filename, TNtuple*& Ntuple, int Debug = 0) {
  size_t buf_size(1000);

  float var[1000];
  int   nvar (-1); // initial

  char  var_names[1000];

  const char delimiters[] = " ,\t\n";
  char* c = (char*) malloc(buf_size);
  
  FILE* f = fopen(Filename,"r");

  Ntuple = NULL;
  
  while (getline(&c,&buf_size,f) >= 0) {
    if (Debug) printf("%s",c);
    if (c[0] == '#') {
//-----------------------------------------------------------------------------
// comment line, could be initialization
//-----------------------------------------------------------------------------
    }
    else {
//-----------------------------------------------------------------------------
// read data, at this point the variable names , if any, should be defined
//-----------------------------------------------------------------------------
      char* token;

      int iv = 0;
      
      token = strtok(c,delimiters);

      parse_number(token,var[iv]);
      if (Debug) printf("token = %s, var = %lf\n",token,var[iv]);
      iv++;
      if (token) {
	while ((token = strtok(NULL,delimiters))) {
	  parse_number(token,var[iv]);
	  if (Debug) printf("token = %s, var = %f\n",token,var[iv]);
	  iv ++;
	}
      }

      if (nvar == -1) {
	nvar = iv;

	var_names[0] = 0;
	for (int i=0; i<nvar; i++) {
	  if (i != 0) strcat(var_names,":");
	  strcat(var_names,Form("var%i",i));
	}
	if (Debug) printf(" INIT: nvar = %i  var_names = \'%s\'\n",nvar,var_names);
//-----------------------------------------------------------------------------
// book ntuple
//-----------------------------------------------------------------------------
	Ntuple = new TNtuple("nt","ntuple",var_names);
      }
      else {
//-----------------------------------------------------------------------------
// make sure each line contains the same number of variables
//-----------------------------------------------------------------------------
	if (iv != nvar) {
	  printf(" ERROR: nvar = %i, iv = %i, BAIL OUT\n",nvar, iv);
	  return -1;
	}
      }

      Ntuple->Fill(var);

      if (Debug) {
	printf(" PARSED INPUT: ");
	for (int i=0; i<nvar; i++) {
	  printf(" %12.5e ",var[i]); 
	}
	printf("\n");
      }
    }
  }
  return 0;
}
