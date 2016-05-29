#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZX.h>
#include <NTL/GF2E.h>
#include <NTL/vec_GF2E.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <list>

using namespace std;
using namespace NTL;

NTL_CLIENT
bool indivisibility(int num_pol, vector<ZZ_pX> pol);
void printTable(vector<ZZ_pX> m, vector<ZZ_pX> b, vector<ZZ_pX> M, vector<ZZ_pX> y, vector<ZZ_pX> bMy);
vector<vector<ZZ_pX>> generateRandomCongruences(int n);
bool isPrime(long long n);
int getMaxDeg(vector<ZZ_pX> polynomials);

int main()
{
	int numOfCongruences;
	int p;

	cout << "Zadajte P pre pole Fp: ";
	cin >> p;

	while (!isPrime(p))
	{
		cout << "Zadali ste zlé číslo. Zadané číslo musí byť prvočíslo." << endl;

		cout << "Zadajte P pre pole Fp: ";
		cin >> p;
	}

	ZZ_p::init(ZZ(p));

	cout << "Koľko kongruencií chcete načítať/vygenerovať (min 2, max 7):  ";
	cin >> numOfCongruences;

	while (numOfCongruences < 2 || numOfCongruences > 7)
    {
		cout << "Zadali ste zlý počet kongruencií." << endl;

		cout << "Koľko kongruencií chcete načítať/vygenerovať (min 2, max 7):  ";
		cin >> numOfCongruences;
	}

	int choice;

	cout << "Chcete zadať kongruencie alebo vygenorovať náhodné?" << endl;
	cout << "1. Zadať" << endl;
	cout << "2. Vygenerovať" << endl;
	cout << "Voľba: ";
	cin >> choice;

	while (choice < 1 || choice > 2)
	{
		cout << "Nesprávna voľba!" << endl;

		cout << "Voľba: ";
		cin >> choice;
	}

	vector<ZZ_pX> m;
	vector<ZZ_pX> b;
	vector<ZZ_pX> M;
	vector<ZZ_pX> y;
	vector<ZZ_pX> bMy;

	if (choice == 1)
	{
		for (int i = 1; i <= numOfCongruences; i++)
		{
			ZZ_pX bi, mi;

			cout << "Zadaj " << i << ". kongruenciu:";
			cin >> bi >> mi;
			m.push_back(mi);
			b.push_back(bi);
		}

		while (!indivisibility(numOfCongruences, m))
		{
			cout << "Polynómi (m) musia byť medzi sebou nesúdeliteľné!" << endl;

			for (int i = 1; i <= numOfCongruences; i++)
			{
				ZZ_pX bi, mi;

				cout << "Zadaj " << i << ". kongruenciu:";
				cin >> bi >> mi;
				m.push_back(mi);
				b.push_back(bi);
			}
		}

		cout << "Vstup O.K." << endl;
	}
	else
	{
		auto congruences = generateRandomCongruences(numOfCongruences);

		m = congruences[0];
		b = congruences[1];
	}

	//VYPOCITA Mi
    for (int i = 0; i < numOfCongruences; i++)
    {
        ZZ_pX Mi;
        bool isMInitialized = false;

        for (int j = 0; j < numOfCongruences; j++)
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
                    Mi *= m[j];
                }
            }
        }

        M.push_back(Mi);
		ZZ_pX invModResult, remainderResult;

		rem(remainderResult, Mi, m[i]);
		InvMod(invModResult, remainderResult, m[i]);

		y.push_back(invModResult);
    }

	for (int i = 0; i < numOfCongruences; i++)
	{
		bMy.push_back(b[i]*M[i]*y[i]);
	}

	//VYPOCITA ZAVERECNE f(x)
	ZZ_pX sum;
	bool sumInitialized = false;
	for (int i = 0; i < numOfCongruences; i++)
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
	for (int i = 0; i < numOfCongruences; i++)
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

	printTable(m, b, M, y, bMy);

	//VYSLEDOK
	ZZ_pX remainderResult;

	rem(remainderResult, sum, MProduct);
	cout << endl << endl << "f(x) = " << remainderResult << endl << endl;

	getchar();
	return 0;
}

vector<vector<ZZ_pX>> generateRandomCongruences(int n)
{
	vector<ZZ_pX> m;
	vector<ZZ_pX> b;

	for (int i = 0; i < n; i++)
	{
		ZZ_pX generated = random_ZZ_pX(8);

		m.push_back(generated);

		int j = 0;

		while (!indivisibility(i + 1, m) && j < 100 || m.back() == conv<ZZ_pX>(1) || m.back() == conv<ZZ_pX>(0))
		{
			m.pop_back();
			generated = random_ZZ_pX(8);
			m.push_back(generated);
			j++;
		}

		if (j == 100)
		{
			m.clear();
			b.clear();
			i = -1;
			continue;
		}

		generated = random_ZZ_pX(8);
		b.push_back(generated);
	}

	auto returnValue = vector<vector<ZZ_pX>>(2);
	returnValue = {m, b};
	return returnValue;
}

//2 2 1 [0 1] [1 1 1] [1 0 1] [1 1 0 1]
bool isPrime(long long n)
{
	long long z = 2;

	while (z * z <= n)
	{
		if (n % z == 0)
			return false;

		z++;
	}

	return true;
}

int getMaxDeg(vector<ZZ_pX> polynomials)
{
	int maxDeg = 0;

	for (int i = 0; i < polynomials.size(); i++)
	{
		auto currentDeg = deg(polynomials[i]);

		if (currentDeg > maxDeg)
		{
			maxDeg = currentDeg;
		}
	}

	return maxDeg;
}

void printTable(vector<ZZ_pX> m, vector<ZZ_pX> b, vector<ZZ_pX> M, vector<ZZ_pX> y, vector<ZZ_pX> bMy)
{
	int columnHeadWidth = 6;
	int columnWidth = getMaxDeg(bMy)*2 + 5;

	cout << endl << endl << left << setw(columnHeadWidth) << "i";
	for (int i = 0; i < m.size(); i++)
	{
		cout << left << setw(columnWidth) << i;
	}

	cout << endl;
	for (int i = 0; i < m.size(); i++)
	{
		cout.fill('-');
		cout << left << setw(columnWidth) << "";
	}
	cout.fill(' ');

	cout << endl << left << setw(columnHeadWidth) << "m";
	for (int i = 0; i < m.size(); i++)
	{
		stringstream ss;
		ss << m[i];
		cout << left << setw(columnWidth) << ss.str();
	}

	cout << endl << left << setw(columnHeadWidth) << "b";
	for (int i = 0; i < b.size(); i++)
	{
		stringstream ss;
		ss << b[i];
		cout << left << setw(columnWidth) << ss.str();
	}

	cout << endl << left << setw(columnHeadWidth) << "M";
	for (int i = 0; i < b.size(); i++)
	{
		stringstream ss;
		ss << M[i];
		cout << left << setw(columnWidth) << ss.str();
	}

	cout << endl << left << setw(columnHeadWidth) << "y";
	for (int i = 0; i < b.size(); i++)
	{
		stringstream ss;
		ss << y[i];
		cout << left << setw(columnWidth) << ss.str();
	}

	cout << endl << left << setw(columnHeadWidth) << "bMy";
	for (int i = 0; i < b.size(); i++)
	{
		stringstream ss;
		ss << bMy[i];
		cout << left << setw(columnWidth) << ss.str();
	}

	cout << endl;
}

bool indivisibility(int num_pol, vector<ZZ_pX> pol)
{
	ZZ_pX result;
	int tmp = num_pol;

	ZZ_pX tmp1, tmp2;

	for (int j = 0; j < num_pol; j++)
    {
        for (int i = 0 + j + 1; i < tmp; i++)
        {
			XGCD(result, tmp1, tmp2, pol[j], pol[i]);
			
			if (conv<ZZ_pX>(1) != result)
				return false;
		}
	}

	return true;
}
