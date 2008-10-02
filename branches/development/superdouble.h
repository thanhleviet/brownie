/*
 * superdouble, a class with double precision but less subject to overflow or underflow
 * superdouble X=mantissa * 10^exponent
 *
 * Brian O'Meara, Oct. 2, 2008
 * GPL2 license
 */
#ifndef SUPERDOUBLE
#define SUPERDOUBLE
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;


class Superdouble  {
private:
	long double mantissa;
	int exponent;
	void adjustDecimal();
	friend ostream& operator<<(ostream& os, const Superdouble& x);
	
public:
		Superdouble(long double mantissa=1.0, int exponent=0);
	~Superdouble();
	Superdouble operator* ( Superdouble x);
	Superdouble operator/ ( Superdouble x);
	Superdouble operator+ ( Superdouble x);
	Superdouble operator- ( Superdouble x);
	void operator++ ();
	void operator -- ();
	void operator*= (Superdouble x);
	void operator/= (Superdouble x);
	void operator+= (Superdouble x);
	void operator-= (Superdouble x);
	int getExponent();
	double getMantissa();
	Superdouble getLn();
//yes, c++ log is log base e, but if log is used here rather than a different word, other calls to log in the program try to use the Superdouble::log rather than the regular log
	
	operator double() {return mantissa*pow(10,exponent);};
	
};


Superdouble::Superdouble(long double m, int e) {
	mantissa=m;
	exponent=e;
	adjustDecimal();
}

Superdouble::~Superdouble() {}

int Superdouble::getExponent() 
{
	return exponent;
}

double Superdouble::getMantissa()
{
	return mantissa;
}

void Superdouble::adjustDecimal() 
{
	if (mantissa==0 || isinf(mantissa) || isnan(mantissa)) {
		exponent=0;
	}
	else {
		while (fabs(mantissa)>=10) {
			mantissa*=0.1;
			exponent+=1;
		}
		while (fabs(mantissa)<1) {
			mantissa*=10.0;
			exponent+=-1;
		}
	}
}

ostream& operator<<(ostream& os, const Superdouble& x)
{
	os<<x.mantissa <<"e"<<x.exponent;
	return os;
}

Superdouble Superdouble::operator * ( Superdouble  x)
{
	Superdouble result(mantissa*x.mantissa,exponent+x.exponent);
	result.adjustDecimal();
	return result;
}

Superdouble Superdouble::operator / ( Superdouble  x)
{
	Superdouble result(mantissa/x.mantissa,exponent-x.exponent);
	result.adjustDecimal();
	return result;
}

Superdouble Superdouble::operator + ( Superdouble  x)
{
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		Superdouble result(mantissa+(x.mantissa*(pow(10,exponentdif))),exponent);
		result.adjustDecimal();
		return result;
	}
	else {
		Superdouble result(x.mantissa,x.exponent);
		result.adjustDecimal();
		return result;
	}
}

Superdouble Superdouble::operator - ( Superdouble  x)
{
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		Superdouble result(mantissa-(x.mantissa*(pow(10,exponentdif))),exponent);
		result.adjustDecimal();
		return result;
	}
	else {
		Superdouble result(-1.0*x.mantissa,x.exponent);
		result.adjustDecimal();
		return result;
	}
}



void Superdouble::operator ++ ()
{
	mantissa++;
	adjustDecimal();
}

void Superdouble::operator -- ()
{
	mantissa--;
	adjustDecimal();
}

void Superdouble::operator *= (Superdouble x)
{
	mantissa*=x.mantissa;
	exponent+=x.exponent;
	adjustDecimal();
}

void Superdouble::operator /= (Superdouble x)
{
	mantissa/=x.mantissa;
	exponent-=x.exponent;
	adjustDecimal();
}

void Superdouble::operator += ( Superdouble  x)
{
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		mantissa=mantissa+(x.mantissa*(pow(10,exponentdif)));
		adjustDecimal();
	}
	else {
		mantissa=x.mantissa;
		exponent=x.exponent;
		adjustDecimal();
	}
}

void Superdouble::operator -= ( Superdouble  x)
{
	//only tricky thing is converting them to same exponent
	if (mantissa!=0) {
		int exponentdif=x.exponent-exponent;
		mantissa=mantissa-(x.mantissa*(pow(10,exponentdif)));
		adjustDecimal();
	}
	else {
		mantissa=-1.0*x.mantissa;
		exponent=x.exponent;
		adjustDecimal();
	}
}


Superdouble Superdouble::getLn () 
{
	//ln(a * 10^b) = ln(a) + ln(10^b) = ln(a) + log10 (10^b) / log10 (e^1) = ln(a) + b/log10(e^1)
	Superdouble result(log(mantissa)+(1.0*(exponent))/log10(exp(1)),0);
	result.adjustDecimal();
	return result;
}

#endif
