#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/tools.h>
#pragma comment (lib,"NTL.lib")
using namespace std;
using namespace NTL;
enum RESULT{MUST_NOT,NOT_SURE,MUST};
ZZ MULTIPLICATIVE_ORDER(ZZ n, ZZ k)/*Ord n(k)*/{
	ZZ ans;
	if (GCD(k, n) != 1){ ans = -1; return ans; }
	else{
		ZZ i;
		for (i = 1;; i++){
			if (PowerMod(k%n, i, n) == 1){ return i; }
		}
	}
}
ZZ EULER(ZZ n){
	ZZ c, i; c = 0;
	for (i = 1; i <= n; i++){
		if (GCD(i, n) == 1)c++;
	}
	return c;
}
RESULT STEP_1(ZZ n){
	ZZ k, s, x0, temp; s = int(log(n) / log(10)) + 1;
	for (k = 2; k <= s; k++){
		x0 = power(to_ZZ("10"), to_long((s - 1) / k + 1));
		while (power(x0, to_long(k)) > n){
			temp = power(x0, to_long(k - 1));
			x0 = x0 - (temp * x0 - n) / (k*temp) - to_long(1);
		}
		if (power(x0, to_long(k)) == n)return MUST_NOT;
	}
	return NOT_SURE;
}
ZZ STEP_2(ZZ n){
	ZZ c; c = sqr(to_ZZ(log(n)/log(2)));//The presicion should be looked after.
	ZZ r;
	for (r = 2;; r++)
		if (MULTIPLICATIVE_ORDER(r, n) > c)return r;
}
RESULT STEP_3(ZZ r, ZZ n){
	ZZ a;
	for (a = 1; a <= r; a++)
		if (GCD(a, n)<n && GCD(a, n)>1)return MUST_NOT;
	return NOT_SURE;
}
RESULT STEP_4(ZZ r, ZZ n){
	if (n<=r)return MUST;
	else return NOT_SURE;
}
bool STEP_5_KERNEL(ZZ r,ZZ a,ZZ n){
	ZZ_p::init(n);
	ZZ_pX m = ZZ_pX(to_long(r), 1) - 1;//m =X^r-1
	ZZ_pX p = ZZ_pX(1, 1) + to_ZZ_p(a);//p =X+a
	ZZ_pX q = ZZ_pX(to_long(n), 1) + to_ZZ_p(a);////q=X^n+a
	ZZ_pX u = PowerMod(p%m, n, m);
	ZZ_pX v = PowerMod(q%m,1,m);
	if (u != v)return false;
	else return true;
}
RESULT STEP_5(ZZ r,ZZ n){
	ZZ c, a; SqrRoot(c,(EULER(r))*log(n));
	for (a = 1; a <= c; a++){
		if (!STEP_5_KERNEL(r, a, n))return MUST_NOT;
	}
	return MUST;
}

int main(){
	ZZ n, r; n = 1; while (cin>>n){
	if (STEP_1(n) == MUST_NOT){cout << "Not prime    "; goto ext; }
	r = STEP_2(n); 
	if (STEP_3(r, n) == MUST_NOT){ cout << "Not prime    "; goto ext; }
	if (STEP_4(r, n) == MUST)    { cout << "Prime        "; goto ext; }
	if (STEP_5(r, n) == MUST)    { cout << "Prime        "; goto ext; }
	else cout << "Not prime    "; ext:; cout << endl;
	}
return 0;
}

