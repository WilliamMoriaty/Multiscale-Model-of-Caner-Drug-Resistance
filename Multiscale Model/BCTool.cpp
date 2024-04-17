#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "time.h"
#include "iostream"

using namespace std;

#include "System.h"
#include "BCTool.h"
#include "Random.h"

struct IMD _MD;
CRandom Rand;

void ReadIPF(char *fn);
void SetParVal(char *str, char const *conststr, char val[]);
void help();
void OutputParameter();

int main(int argc, char *argv[])
{
	System *sys;
	if (argc < 2)
	{
		help();
		exit(0);
	}

	ReadIPF(argv[1]);
	// OutputParameter();
	if (_MD.seed > 0)
		Rand.Initialized(_MD.seed);
	sys = new System(_MD.N0);
	if (sys->Initialized())
	{
		sys->Run();
		// OutputParameter();
		return (1);
	}
	else
	{
		return (0);
	}
}
void help()
{
	cout << "Usage: " << endl;
	cout << "bct_VirusCell inputfile" << endl;
	cout << "inputfile: The name of input file." << endl;
	cout << " For example: bct md.in." << endl;
}
void SetParVal(char *str, char const *conststr, char val[])
{
	char *pst;
	if ((pst = strstr(str, conststr)) != NULL)
	{
		strcpy(val, pst + strlen(conststr));
		if ((pst = strstr(val, "\"")) != NULL)
		{
			val[pst - val] = '\0';
		}
	}
	return;
}
void ReadIPF(char *fn)
{
	FILE *fp, *fp1;
	char str[Strlength], *pst;
	if ((fp = fopen(fn, "r")) == NULL)
	{
		cout << "Cannot open the input file." << endl;
		exit(0);
	}
	rewind(fp);
	while (!feof(fp))
	{
		fgets(str, Strlength, fp);
		if (str[0] == '#')
		{
			continue;
		}
		SetParVal(str, "mdcrd=\"", _MD.mdcrd);
		SetParVal(str, "cellpar=\"", _MD.cellpar);
		SetParVal(str, "dose=\"", _MD.dose);
		if ((pst = strstr(str, "dt=")) != NULL)
		{
			_MD.dt = atof(pst + 3);
		}
		if ((pst = strstr(str, "T1=")) != NULL)
		{
			_MD.T1 = atof(pst + 3);
		}
		if ((pst = strstr(str, "T0=")) != NULL)
		{
			_MD.T0 = atof(pst + 3);
		}
		if ((pst = strstr(str, "ntpx=")) != NULL)
		{
			_MD.ntpx = atoi(pst + 5);
		}
		if ((pst = strstr(str, "ntpr=")) != NULL)
		{
			_MD.ntpr = atoi(pst + 5);
		}
		if ((pst = strstr(str, "seed=")) != NULL)
		{
			_MD.seed = atoi(pst + 5);
		}
		if ((pst = strstr(str, "N0=")) != NULL)
		{
			_MD.N0 = atoi(pst + 3);
		}
	}
	// load dose c
	if ((fp1 = fopen(_MD.dose, "r")) == NULL)
	{
		cout << "Cannot open the input dose file." << endl;
		exit(0);
	}
	for (int i = 0; i < 7200; i++)
	{
		fscanf(fp1, "%f\n", &_MD.c[i]);
	}

	fclose(fp);
	fclose(fp1);
	return;
}
void OutputParameter()
{

	cout << "x=" << _MD.c << endl;
	
	cout << "T1=" << _MD.T1 << endl;
	//	 cout<<"ntpr="<<_MD.ntpr<<endl;
	//	 cout<<"ntpx="<<_MD.ntpx<<endl;
}
