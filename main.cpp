// NTLtest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZX.h>
#include <cassert>
#include "NTL/ZZ.h"
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>
#include <NTL/ZZX.h>
using namespace std;
using namespace NTL;

NTL_CLIENT
void GCDx(ZZX& d, const ZZX& a, const ZZX& b);
bool indivisibility(int num_pol, ZZ_pX *pol);

int main()
{
	
	Vec<ZZ_p> vec_ZZ_p;
	//ZZ_p::init(ZZ(3)); // define GF(p)
	//ZZ_pX z;

	int num_pol;	//number of polynomials
	int f;			//field or ring
	cout << "Definuj P:";
	cin >> f;
	ZZ_p::init(ZZ(f)); // define GF(p)
				//polynom

	cout << "[min 2 , max 7]Kolko polynomov chces nacitat: ";
	cin >> num_pol;
	
	ZZ_pX *pol=new ZZ_pX[num_pol];			//array of polynomials

	if (num_pol < 1){
		cout << "Zadali ste zly vstup to je malo";
		return 0;
	}

	if (num_pol > 7){
		cout << "Zadali ste zly vstup to je vela";
		return 0;
	}
	
	for (int i = 1; i <= num_pol; i++){
		cout << "Zadaj " << i << ". polynom:";
		cin >> pol[i-1];
	}
	
	if (indivisibility(num_pol, pol))         // function verfies indivisibility of polynomials in array
		cout << "Su delitelne medzi sebou";
	else
		cout << "Niesu delitelne medzi sebou";



	getchar();
	return 0;
}


bool indivisibility(int num_pol, ZZ_pX *pol){
	ZZX result;
	double vOne = 1;		//value 1
	int tmp = num_pol;

	for (int j = 0; j < num_pol; j++){

		for (int i = 0 + j; i < tmp; i++){

			GCDx(result, conv<ZZX>(pol[j]), conv<ZZX>(pol[i]));	//gcd of two polynomials
			
			if (conv<ZZX>(vOne) == result)			//test if GCD of polynomials are 1
				return true;
		}
	}

	return false;
}



void GCDx(ZZX& d, const ZZX& a, const ZZX& b)
{
	if (a == 0) {
		d = b;
		if (LeadCoeff(d) < 0) negate(d, d);
		return;
	}

	if (b == 0) {
		d = a;
		if (LeadCoeff(d) < 0) negate(d, d);
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
	for (i = 0;; i++) {
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


		if (deg(G) == 0) {
			res = 1;
			break;
		}

		if (FirstTime || deg(G) < deg(g)) {
			prod = 1;
			g = 0;
			FirstTime = 0;
		}
		else if (deg(G) > deg(g)) {
			continue;
		}

		if (!CRT(g, prod, G)) {
			PrimitivePart(res, g);
			if (divide(f1, res) && divide(f2, res))
				break;
		}

	}

	mul(d, res, c);
	if (LeadCoeff(d) < 0) negate(d, d);
}
