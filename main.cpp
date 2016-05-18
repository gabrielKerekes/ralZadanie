#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZX.h>
#include <cassert>
#include "NTL/ZZ.h"
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>
#include <NTL/ZZX.h>
#include <vector>

using namespace std;
using namespace NTL;

NTL_CLIENT
void GCDx(ZZX& d, const ZZX& a, const ZZX& b);
bool indivisibility(int num_pol, vector<ZZ_pX> pol);

int main()
{
	Vec<ZZ_p> vec_ZZ_p;
	//ZZ_p::init(ZZ(3)); // define GF(p)
	//ZZ_pX z;

	int num_pol;	//number of polynomials
	int f;			//field or ring
	cout << "Definuj P:";
	cin >> f;
	ZZ_p::init(ZZ(f));

	cout << "[min 2 , max 7] Kolko kongruencii chces nacitat: ";
	cin >> num_pol;
	
	ZZ_pX *pol = new ZZ_pX[num_pol];			//array of polynomials
    vector<ZZ_pX> m;
    vector<ZZ_pX> b;
    vector<ZZ_pX> M;
    vector<ZZ_pX> y;
    vector<ZZ_pX> bMy;

	if (num_pol < 1)
    {
		cout << "Zadali ste zly vstup to je malo";
		return 0;
	}

	if (num_pol > 7)
    {
		cout << "Zadali ste zly vstup to je vela";
		return 0;
	}
	
	for (int i = 1; i <= num_pol; i++)
    {
        ZZ_pX bi, mi;

		cout << "Zadaj " << i << ". kongruenciu:";
		//cin >> pol[i-1];
        cin >> bi >> mi;
        m.push_back(mi);
        b.push_back(bi);
	}
	
	if (indivisibility(num_pol, m))         // function verfies indivisibility of polynomials in array
		cout << "Su delitelne medzi sebou" << endl;
	else
		cout << "Niesu delitelne medzi sebou" << endl;

	//VYPOCITA Mi
    for (int i = 0; i < num_pol; i++)
    {
        ZZ_pX Mi;
        bool isMInitialized = false;

        for (int j = 0; j < num_pol; j++)
        {
            if (i != j)
            {
                if (!isMInitialized)
                {
                    Mi = m[j];
                    isMInitialized = true;
                }
                else
                {
                    Mi *= M[j];
                }
            }
        }

        M.push_back(Mi);
		//cout << Mi << " " << m[i] << endl;
		ZZ_pX invModResult, remainderResult;
		//cout << endl << "DEG M: " << deg(Mi) << " DEG m[i]: " << deg(m[i]) << endl;
		rem(remainderResult, Mi, m[i]);
		//cout << "REM: " << remainderResult << endl;
		InvMod(invModResult, remainderResult, m[i]);
		//cout << endl << "INVMOD: " << invModResult << endl;
		y.push_back(invModResult);
    }

	cout << "M: ";
    for (int i = 0; i < num_pol; i++)
    {
        cout << M[i] << " ";
    }
	cout << endl;

	cout << "y: ";
	for (int i = 0; i < num_pol; i++)
	{
		cout << y[i] << " ";
	}
	cout << endl;

	for (int i = 0; i < num_pol; i++)
	{
		bMy.push_back(b[i]*M[i]*y[i]);
	}

	cout << "y: ";
	for (int i = 0; i < num_pol; i++)
	{
		cout << bMy[i] << " ";
	}
	cout << endl;

	//VYPOCITA ZAVERECNE f(x)
	ZZ_pX sum;
	bool sumInitialized = false;
	for (int i = 0; i < num_pol; i++)
	{
		if (!sumInitialized)
		{
			sum = bMy[i];
			sumInitialized = true;
		}
		else
		{
			sum += bMy[i];
		}
	}
	//VYPOCITA M
	ZZ_pX MProduct;
	bool MProductInitialized = false;
	for (int i = 0; i < num_pol; i++)
	{
		if (!MProductInitialized)
		{
			MProduct = m[i];
			MProductInitialized = true;
		}
		else
		{
			MProduct *= m[i];
		}
	}
	//VYSLEDOK
	ZZ_pX remainderResult;
	//cout << endl << "DEG M: " << deg(Mi) << " DEG m[i]: " << deg(m[i]) << endl;
	rem(remainderResult, sum, MProduct);
	cout << endl << endl << "FINAL: " << remainderResult << endl << endl;

	getchar();
	return 0;
}

bool indivisibility(int num_pol, vector<ZZ_pX> pol)
{
	ZZX result;
	double vOne = 1;		//value 1
	int tmp = num_pol;

	for (int j = 0; j < num_pol; j++)
    {
        for (int i = 0 + j; i < tmp; i++)
        {
			GCDx(result, conv<ZZX>(pol[j]), conv<ZZX>(pol[i]));	//gcd of two polynomials
			
			if (conv<ZZX>(vOne) == result)			//test if GCD of polynomials are 1
				return true;
		}
	}

	return false;
}

void GCDx(ZZX& d, const ZZX& a, const ZZX& b)
{
	if (a == 0)
    {
		d = b;
		if (LeadCoeff(d) < 0) NTL::negate(d, d);
		return;
	}

	if (b == 0)
    {
		d = a;
		if (LeadCoeff(d) < 0) NTL::negate(d, d);
		return;
	}

	ZZ c1, c2, c;
	ZZX f1, f2;

	content(c1, a);
	divide(f1, a, c1);

	content(c2, b);
	divide(f2, b, c2);

	GCD(c, c1, c2);

	ZZ ld;
	GCD(ld, LeadCoeff(f1), LeadCoeff(f2));

	ZZX g, res;

	ZZ prod;

	zz_pPush push; // save current modulus, restore upon return

	long FirstTime = 1;

	long i;
	for (i = 0;; i++)
    {
		zz_p::FFTInit(i);
		long p = zz_p::modulus();

		if (divide(LeadCoeff(f1), p) || divide(LeadCoeff(f2), p)) continue;

		zz_pX G, F1, F2;
		zz_p  LD;

		conv(F1, f1);
		conv(F2, f2);
		conv(LD, ld);

		GCD(G, F1, F2);
		mul(G, G, LD);

		if (deg(G) == 0)
        {
			res = 1;
			break;
		}

		if (FirstTime || deg(G) < deg(g))
        {
			prod = 1;
			g = 0;
			FirstTime = 0;
		}
		else if (deg(G) > deg(g))
        {
			continue;
		}

		if (!CRT(g, prod, G))
        {
			PrimitivePart(res, g);
			if (divide(f1, res) && divide(f2, res))
				break;
		}
	}

	mul(d, res, c);
	if (LeadCoeff(d) < 0) NTL::negate(d, d);
}
