#define INI_MAX 1000
char * ini_key[INI_MAX];
char * ini_val[INI_MAX];
int ini_total = 0;
int nReadIniNumber=0;
int nGeomTransfer=0;
#include <cstring>
#include <stdarg.h>

//return a string without leading and tailing spaces
char * trim(char * buf)
{
	int iStart = 0, iEnd = strlen(buf) - 1;

	for (; iStart <= iEnd; iStart++)
	{
		char s = buf[iStart];
		if (s != ' ' && s != '\n' && s != '\r') break;
	}

	//skip trailing spaces
	for (; iStart <= iEnd; iEnd--)
	{
		char s = buf[iEnd];
		if (s != ' ' && s != '\n' && s != '\r') break;
	}

	buf[iEnd + 1] = 0;
	return buf + iStart;
}

//log the error to "error.log"
void errorLog(const char * sFormat, ...)
{
	FILE * fp;
	time_t ltime;
	va_list lst;
	va_start(lst, sFormat);

	printf("[error]");
	vprintf(sFormat, lst);
	printf("\n");

	if ((fp = fopen("error.log", "a")) != NULL)
	{
		time(&ltime);
		fprintf(fp, "\n%s", ctime(&ltime));

		vfprintf(fp, sFormat, lst);
		fclose(fp);
	}
}


int IniParse(const char * fname)
{
	ini_total=0;
	FILE * fp;
	char buf[256], *cpos, *sKey, *sVal;
	char sGroup[256] = { "" };

	if ((fp = fopen(fname, "r")) != NULL)
	{
		while (!feof(fp))
		{
			if (fgets(buf, sizeof(buf), fp) == NULL)
			{
				if (ferror(fp))
				{
					sprintf(buf, "ParseIni: failed in %s", fname);
					errorLog(buf);
					return -1;
				}
				break;
			}
			if (buf[0] == '#') continue;	//comment line

			cpos = strchr(buf, '=');
			if (cpos != NULL)
			{	//key=value

				*cpos++ = 0;	//break the string
				sKey = trim(buf);

				sVal = trim(cpos);
				
				if (ini_total < INI_MAX)
				{
					ini_key[ini_total] = new char[strlen(sKey) + 1];
					strcpy(ini_key[ini_total], sKey);
					ini_val[ini_total] = new char[strlen(sVal) + 1];
					strcpy(ini_val[ini_total], sVal);
					ini_total++;
				}
			}
		}
		fclose(fp);
		return 0;
	}

	return -1;
}

void IniUnload()
{
	for (int i = 0; i < ini_total; i++)
	{
		delete ini_key[i];
		delete ini_val[i];
	}
}

char * IniFindValue(const char * key)
{
	for (int i = 0; i < ini_total; i++)
	{
		if( strcmp(key, ini_key[i]) == 0 )
		return ini_val[i];
	}
	return NULL;
}



int ReadIni(const char* fname, double &lx, double &ly, int &Nx, int &Ny, 
			int &variant, double &tmax, double &cu, int &Rodionov, int &tout)
{

	IniParse(fname);
	char *p;

	p = IniFindValue("lx");
	if(p != NULL)
	{
		lx = atof(p);
	}

	p = IniFindValue("ly");
	if(p != NULL)
	{
		ly = atof(p);
	}

	p = IniFindValue("Nx");
	if(p != NULL)
	{
		Nx = atoi(p);
	}

	p = IniFindValue("Ny");
	if(p != NULL)
	{
		Ny = atoi(p);
	}

	p = IniFindValue("variant");
	if(p != NULL)
	{
		variant = atoi(p);
	}

	p = IniFindValue("tmax");
	if(p != NULL)
	{
		tmax = atof(p);
	}

	p = IniFindValue("cu");
	if(p != NULL)
	{
		cu = atof(p);
	}

	p = IniFindValue("Rodionov");
	if(p != NULL)
	{
		Rodionov = atoi(p);
	}

	p = IniFindValue("tout");
	if(p != NULL)
	{
		tout = atoi(p);
	}

	IniUnload();
	return 0;
}