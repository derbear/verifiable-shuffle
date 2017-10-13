/*
 * func_ver.cpp
 *
 *  Created on: 04.07.2012
 *      Author: stephaniebayer
 */

#include "func_ver.h"
#include<vector>
#include "Cipher_elg.h"
#include "G_q.h"
#include "Mod_p.h"
#include "Functions.h"
#include "ElGammal.h"
#include "multi_expo.h"
#include <fstream>

#include <time.h>
#include <NTL/ZZ.h>
NTL_CLIENT


extern G_q G;
extern G_q H;
extern Pedersen Ped;
extern ElGammal El;

func_ver::func_ver() {
	// TODO Auto-generated constructor stub

}

func_ver::~func_ver() {
	// TODO Auto-generated destructor stub
}


int func_ver::check_Dh(vector<Mod_p>* c_Dh, vector<ZZ>* chal, vector<ZZ>* D_h_bar, ZZ r_Dh_bar){
	long i;
	Mod_p t_Dh, co_Dh, temp;
	long m = c_Dh->size();

	t_Dh =c_Dh->at(0);
	for (i=1; i<m;i++){
		Mod_p::expo(temp,c_Dh->at(i),chal->at(i-1));
		Mod_p::mult(t_Dh,t_Dh,temp);
	}
	co_Dh = Ped.commit(D_h_bar,r_Dh_bar);
/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"D_h "<<t_Dh<<" "<<co_Dh<<endl;*/
	if (t_Dh == co_Dh){
		return 1;
	}
	return 0;
}

int func_ver::check_D(Mod_p c_D0, Mod_p c_z, vector<Mod_p>* c_A, vector<Mod_p>* c_B, vector<ZZ>* chal_1, ZZ chal_2, vector<ZZ>* A_bar, ZZ r_A_bar, long n){
	int i;
	Mod_p t_D, co_D, temp, inv;
	long m = c_A->size()-1;
	ZZ ord = H.get_ord();
	vector<ZZ>* v_1 = new vector<ZZ>(n);

	t_D = c_D0;
	Mod_p::inv(inv, c_z);
	for (i=1; i<m;i++){
		Mod_p::expo(temp,c_A->at(i),chal_2);
		Mod_p::mult(temp, temp,c_B->at(i));
		Mod_p::mult(temp,temp,inv);
		Mod_p::expo(temp, temp, chal_1->at(i-1));
		Mod_p::mult(t_D,t_D,temp);
	}
	for(i=0; i<n;i++){
		NegateMod(v_1->at(i),to_ZZ(1),ord);
	}
	temp=Ped.commit(v_1,to_ZZ(0));
	Mod_p::expo(temp, temp, chal_1->at(m-1));
	Mod_p::mult(t_D,t_D,temp);
	co_D = Ped.commit(A_bar, r_A_bar);

	/*string name = "example.txt";
	ofstream s;
	ost.open(name.c_str(),ios::app);
	ost<<"D "<<t_D<<" "<<co_D<<endl;*/
	delete v_1;
	if(t_D == co_D){
		return 1;
	}
	return 0;
}


int func_ver::check_Ds(vector<Mod_p>* c_Ds, vector<Mod_p>* c_Dh, Mod_p c_Dm, vector<ZZ>* chal_1, vector<ZZ>* chal_2, vector<ZZ>* Ds_bar, ZZ r_Ds_bar){
	long i,l;
	long m= c_Ds->size()-1;
	Mod_p t_Ds, co_Ds, temp, temp_1;

	l=m-1;
	for(i=0; i<l; i++){
		Mod_p::expo(c_Ds->at(i),c_Dh->at(i), chal_1->at(i));
	}
	if(m>1){
		Mod_p::expo(temp, c_Dh->at(1), chal_1->at(0));
		for(i=1;i<l; i++){
			Mod_p::expo(temp_1,c_Dh->at(i+1), chal_1->at(i));
			Mod_p::mult(temp,temp,temp_1);
		}
		c_Ds->at(l)=temp;
	}
	else{
		c_Ds->at(l)=Mod_p(1,G.get_mod());
	}

	c_Ds->at(m)=c_Dm;
	Mod_p::expo(t_Ds, c_Ds->at(0),chal_2->at(m-1));
	for(i=1; i<m; i++){
		Mod_p::expo(temp, c_Ds->at(i), chal_2->at(m-1-i));
		Mod_p::mult(t_Ds, t_Ds,temp);
	}
	Mod_p::mult(t_Ds,t_Ds,c_Ds->at(m));
	co_Ds = Ped.commit(Ds_bar, r_Ds_bar);

/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"Ds "<<t_Ds<<" "<<co_Ds<<endl;*/
	if(t_Ds == co_Ds){
		return 1;
	}
	return 0;
}

int func_ver::check_Dl(vector<Mod_p>* c_Dl, vector<ZZ>* chal_1, vector<ZZ>* A, vector<ZZ>* B, vector<ZZ>*  chal_2, ZZ r_Dl_bar){
	long i;
	Mod_p t_Dl, co_Dl, temp;
	ZZ temp_1;
	long l = c_Dl->size();
	long pos = (l-1)/2+1;
	ZZ mod = G.get_mod();
	t_Dl = c_Dl->at(0);
	for(i=1; i<l; i++){
		Mod_p::expo(temp, c_Dl->at(i),chal_1->at(i-1));
		Mod_p::mult(t_Dl, t_Dl, temp);
	}
	temp_1=Functions::bilinearMap(A,B,chal_2);
	co_Dl =Ped.commit(temp_1,r_Dl_bar);

/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"Dl "<<t_Dl<<" "<<co_Dl<<endl;*/

	temp = Mod_p(1,mod);
	if(t_Dl==co_Dl & c_Dl->at(pos) ==temp){
		return 1;
	}
	return 0;
}


int func_ver::check_d(vector<Mod_p>* c_Dh, Mod_p c_d, vector<ZZ>* chal, vector<ZZ>* d_bar, ZZ r_d_bar){
	Mod_p t_d, co_d, temp;
	long m = c_Dh->size();

	Mod_p::expo(temp, c_Dh->at(m-1), chal->at(0));
	Mod_p::mult(t_d, temp, c_d);
		co_d = Ped.commit(d_bar, r_d_bar);

	/*	string name = "example.txt";
		ofstream ost;
		ost.open(name.c_str(),ios::app);
		ost<<"d "<<t_d<<" "<<co_d<<endl;*/
	if(t_d==co_d){
		return 1;
	}
	return 0;
}


int func_ver::check_Delta(Mod_p c_dh, Mod_p c_Delta, vector<ZZ>* chal, vector<ZZ>* Delta_bar, vector<ZZ>* d_bar, ZZ r_Delta_bar, ZZ chal_1, ZZ chal_2, ZZ chal_3){
	long i,j;
	Mod_p t_Delta, co_Delta, temp;
	ZZ t_1, t_2, t_3, prod, chal_temp;
	ZZ ord = H.get_ord();
	long m = (chal->size()-1)/2;
	long n = Delta_bar->size();
	vector<ZZ>* Delta_temp=0;

	Mod_p::expo(temp, c_dh, chal->at(0));
	Mod_p::mult(t_Delta, temp, c_Delta);

	Delta_temp = new vector<ZZ>(n-1);
	t_3= chal->at(0);
	for(i=0; i<n-1; i++){
		MulMod(t_1, Delta_bar->at(i), d_bar->at(i+1), ord);
		MulMod(t_2, t_3, Delta_bar->at(i+1),ord);
		SubMod(Delta_temp->at(i), t_2, t_1, ord);
	}

	co_Delta = Ped.commit(Delta_temp, r_Delta_bar);

	delete Delta_temp;
/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"Delta "<<t_Delta<<" "<<c_Delta<<endl;*/

	prod = to_ZZ(1);
	chal_temp =to_ZZ(1);
	for(i=1; i<=m; i++){
		for(j=1; j<=n; j++){
			MulMod(chal_temp, chal_temp, chal_1,ord);
			SubMod(t_1, chal_temp, chal_2,ord);
			t_3 = n*(i-1)+j;
			MulMod(t_3,t_3, chal_3, ord);
			AddMod(t_1,t_1, t_3, ord);
			MulMod(prod,prod, t_1, ord);
		}
	}
	MulMod(prod, prod, chal->at(0), ord);

	//	ost<<"prod "<<prod<<" "<<Delta_bar->at(n-1)<<endl;
	//ost<<d_bar->at(0)<<" "<<Delta_bar->at(0)<<endl;
	if(t_Delta ==co_Delta)
		if(prod ==Delta_bar->at(n-1) & d_bar->at(0)==Delta_bar->at(0)){
		return 1;
	}
	return 0;
}

int func_ver::check_B(vector<Mod_p>* c_B, Mod_p c_B0, vector<ZZ>* chal, vector<ZZ>* B_bar, ZZ r_B_bar){
	long i;
	Mod_p t_B, co_B, temp, temp_1;
	long m = c_B->size();
	//check for correctness of the committed B
	temp = c_B0;
	for(i = 0; i<m; i++){
		Mod_p::expo(temp_1,c_B->at(i), chal->at(i));
		Mod_p::mult(temp,temp,temp_1);
	}
	t_B = temp;
	co_B = Ped.commit(B_bar, r_B_bar);
/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"B "<<t_B<<" "<<co_B<<endl;*/
	if(t_B==co_B){
		return 1;
	}
	return 0;
}

int func_ver::check_a(vector<Mod_p>* c_a, vector<Mod_p>* c_Dl, vector<ZZ>* chal, ZZ a_bar, ZZ r_a_bar){
	long i,l;
	Mod_p t_a, co_a, temp, temp_1;
	long m = c_a->size()/2;

	//Check that the random values are used right
	temp = c_a->at(0);
	l=2*m-1;
	for(i = 1; i<=l; i++){
		Mod_p::expo(temp_1,c_a->at(i), chal->at(i-1));
		Mod_p::mult(temp,temp,temp_1);
	}
	t_a = temp;
	co_a = Ped.commit(a_bar, r_a_bar);
/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"a "<<t_a<<" "<<co_a<<endl;
*/
	if(t_a==co_a & c_a->at(m)==c_Dl->at(m+1)){//both commitments should be com(0,0), so equal
		return 1;
	}
	return 0;
}

int func_ver::check_c(vector<vector<Cipher_elg>* >* enc, vector<Cipher_elg>* E, ZZ chal){
	long i,j;
	ZZ chal_temp;
	ZZ ord = H.get_ord();
	Cipher_elg c, temp;
	long m = enc->size();
	long n = enc->at(0)->size();

	chal_temp = to_ZZ(1);
	c = Cipher_elg(1,1,H.get_mod());
 	for(i = 0; i<m;i++){
		for(j = 0; j<n;j++){
			MulMod(chal_temp, chal_temp, chal, ord);
			Cipher_elg::expo(temp,enc->at(i)->at(j),chal_temp);
			Cipher_elg::mult(c,c,temp);
		}
	}
 	if(E->at(m)==c){
 		return 1;
 	}
 	return 0;
}

int func_ver::check_E(vector<vector<Cipher_elg>* >* C, vector<Cipher_elg>* E, vector<ZZ>* chal, vector<ZZ>* B_bar, ZZ a_bar, ZZ rho_bar ){
	long i,j,l;
	Cipher_elg t_E, co_E, temp, temp_1, temp_2;
	Mod_p gen,t;
	long m = C->size();
	long n = C->at(0)->size();
	ZZ te;
	ZZ ord = H.get_ord();
	ZZ mod = H.get_mod();

 	temp = E->at(0);
 	l=2*m;
	for(i = 1; i<l; i++){
		Cipher_elg::expo(temp_1,E->at(i), chal->at(i-1));
		Cipher_elg::mult(temp,temp,temp_1);
	}
	t_E = temp;

	gen = H.get_gen();
	t = Mod_p::expo(gen,a_bar);
	temp = El.encrypt(t, rho_bar);
	temp_1=Cipher_elg(1,1,mod);
	l=m-1;
	for(i = 0; i<l;i++){
		for(j = 0; j<n;j++){
			MulMod(te , B_bar->at(j),chal->at(m-i-2),ord);
			Cipher_elg::expo(temp_2,C->at(i)->at(j),te);
			Cipher_elg::mult(temp_1,temp_1,temp_2);
		}
	}
	for(j = 0; j<n;j++){
		te =B_bar->at(j);
		Cipher_elg::expo(temp_2,C->at(m-1)->at(j),te);
		Cipher_elg::mult(temp_1,temp_1,temp_2);
	}
	Cipher_elg::mult(co_E ,temp_1, temp);

/*	string name = "example.txt";
	ofstream ost;
	ost.open(name.c_str(),ios::app);
	ost<<"E"<<t_E<<" "<<co_E<<endl;*/
	if(t_E==co_E){
		return 1;
	}
	return 0;
}


int func_ver::check_Dh_op(vector<Mod_p>* c_Dh, vector<ZZ>* chal, vector<ZZ>* D_h_bar, ZZ r_Dh_bar, long win_LL){
	Mod_p t_Dh, co_Dh;

	multi_expo::multi_expo_LL(t_Dh,c_Dh, chal, win_LL);
	co_Dh = Ped.commit_opt(D_h_bar,r_Dh_bar);

	//cout<<"D_h "<<t_Dh<<" "<<co_Dh<<endl;
	if (t_Dh == co_Dh){
		return 1;
	}
	return 0;
}


int func_ver::check_D_op(Mod_p c_D0, Mod_p c_z, vector<Mod_p>* c_A, vector<Mod_p>* c_B, vector<ZZ>* chal_1, ZZ chal_2, vector<ZZ>* A_bar, ZZ r_A_bar, long n){
	int i;
	Mod_p t_D, co_D, temp, inv;
	long m = c_A->size()-1;
	ZZ ord = H.get_ord();
	vector<ZZ>* v_1 = new vector<ZZ>(n);

	t_D = c_D0;
	Mod_p::inv(inv, c_z);
	for (i=1; i<m;i++){
		Mod_p::expo(temp,c_A->at(i),chal_2);
		Mod_p::mult(temp, temp,c_B->at(i));
		Mod_p::mult(temp,temp,inv);
		Mod_p::expo(temp, temp, chal_1->at(i-1));
		Mod_p::mult(t_D,t_D,temp);
	}
	for(i=0; i<n;i++){
		NegateMod(v_1->at(i),to_ZZ(1),ord);
	}
	temp=Ped.commit_opt(v_1,to_ZZ(0));
	Mod_p::expo(temp, temp, chal_1->at(m-1));
	Mod_p::mult(t_D,t_D,temp);
	co_D = Ped.commit_opt(A_bar, r_A_bar);
	//cout<<"D "<<t_D<<" "<<co_D<<endl;
	delete v_1;
	if(t_D == co_D){
		return 1;
	}
	return 0;
}


int func_ver::check_Ds_op(vector<Mod_p>* c_Ds, vector<Mod_p>* c_Dh, Mod_p c_Dm, vector<ZZ>* chal_1, vector<ZZ>* chal_2, vector<ZZ>* Ds_bar, ZZ r_Ds_bar){
	long i,l;
	long m= c_Ds->size()-1;
	Mod_p t_Ds, co_Ds, temp, temp_1;

	l=m-1;
	for(i=0; i<l; i++){
		Mod_p::expo(c_Ds->at(i),c_Dh->at(i), chal_1->at(i));
	}
	if(m>1){
		Mod_p::expo(temp, c_Dh->at(1), chal_1->at(0));
		for(i=1;i<l; i++){
			Mod_p::expo(temp_1,c_Dh->at(i+1), chal_1->at(i));
			Mod_p::mult(temp,temp,temp_1);
		}
		c_Ds->at(l)=temp;
	}
	else{
		c_Ds->at(l)=Mod_p(1,G.get_mod());
	}

	c_Ds->at(m)=c_Dm;
	Mod_p::expo(t_Ds, c_Ds->at(0),chal_2->at(m-1));
	for(i=1; i<m; i++){
		Mod_p::expo(temp, c_Ds->at(i), chal_2->at(m-1-i));
		Mod_p::mult(t_Ds, t_Ds,temp);
	}
	Mod_p::mult(t_Ds,t_Ds,c_Ds->at(m));
	co_Ds = Ped.commit_opt(Ds_bar, r_Ds_bar);
	//cout<<"Ds "<<t_Ds<<" "<<co_Ds<<endl;
	if(t_Ds == co_Ds){
		return 1;
	}
	return 0;
}

int func_ver::check_Dl_op(vector<Mod_p>* c_Dl, vector<ZZ>* chal, vector<ZZ>* A_bar, vector<ZZ>* Ds_bar, vector<ZZ>*  chal_1, ZZ r_Dl_bar){
	long i;
	Mod_p t_Dl, co_Dl, temp;
	ZZ temp_1;
	long l = c_Dl->size();
	long pos = (l-1)/2+1;
	ZZ mod = G.get_mod();

	t_Dl = c_Dl->at(0);
	for(i=1; i<l; i++){
		Mod_p::expo(temp, c_Dl->at(i),chal->at(i-1));
		//Mod_p::mult(t_Dl,t_Dl, temp);
		t_Dl = t_Dl*temp;
	}
	temp_1=Functions::bilinearMap(A_bar,Ds_bar,chal_1);
	co_Dl =Ped.commit_sw(temp_1,r_Dl_bar);
	//cout<<"Dl "<<t_Dl<<" "<<co_Dl<<endl;

	temp= Mod_p(1, mod);
	if(t_Dl==co_Dl & c_Dl->at(pos)==temp){
		return 1;
	}
	return 0;
}


int func_ver::check_d_op(vector<Mod_p>* c_Dh, Mod_p c_d, vector<ZZ>* chal, vector<ZZ>* d_bar, ZZ r_d_bar){
	Mod_p t_d, co_d, temp;
	long m = c_Dh->size();
	Mod_p::expo(temp, c_Dh->at(m-1), chal->at(0));
	Mod_p::mult(t_d, temp, c_d);
	co_d = Ped.commit_opt(d_bar, r_d_bar);
	//cout<<"d "<<t_d<<" "<<co_d<<endl;
	if(t_d==co_d){
		return 1;
	}
	return 0;
}


int func_ver::check_Delta_op(Mod_p c_dh, Mod_p c_Delta, vector<ZZ>* chal, vector<ZZ>* Delta_bar, vector<ZZ>* d_bar, ZZ r_Delta_bar, ZZ chal_1, ZZ chal_2, ZZ chal_3){
	long i,j;
	Mod_p t_Delta, co_Delta, temp;
	ZZ t_1, t_2, t_3, prod, chal_temp;
	ZZ ord = H.get_ord();
	long m = (chal->size()-1)/2;
	long n = Delta_bar->size();
	vector<ZZ>* Delta_temp=0;

	Mod_p::expo(temp, c_dh, chal->at(0));
	Mod_p::mult(t_Delta, temp, c_Delta);

	Delta_temp = new vector<ZZ>(n-1);
	t_3= chal->at(0);
	for(i=0; i<n-1; i++){
		MulMod(t_1, Delta_bar->at(i), d_bar->at(i+1), ord);
		MulMod(t_2, t_3, Delta_bar->at(i+1),ord);
		SubMod(Delta_temp->at(i), t_2, t_1, ord);
	}

	co_Delta = Ped.commit_opt(Delta_temp, r_Delta_bar);

	delete Delta_temp;
	//cout<<"Delta "<<t_Delta<<" "<<co_Delta<<endl;

	prod = to_ZZ(1);
	chal_temp =to_ZZ(1);
	for(i=1; i<=m; i++){
		for(j=1; j<=n; j++){
			MulMod(chal_temp, chal_temp, chal_1,ord);
			SubMod(t_1, chal_temp, chal_2,ord);
			t_3 = n*(i-1)+j;
			MulMod(t_3,t_3, chal_3, ord);
			AddMod(t_1,t_1, t_3, ord);
			MulMod(prod,prod, t_1, ord);
		}
	}
	MulMod(prod, prod, chal->at(0), ord);

	//cout<<"prod "<<prod<<" "<<Delta_bar->at(n-1)<<endl;
	//cout<<d_bar->at(0)<<" "<<Delta_bar->at(0)<<endl;
	if(t_Delta ==co_Delta)
		if(prod ==Delta_bar->at(n-1) & d_bar->at(0)==Delta_bar->at(0)){
		return 1;
	}
	return 0;
}

int func_ver::check_B_op(vector<Mod_p>* c_B, Mod_p c_B0, vector<ZZ>* chal, vector<ZZ>* B_bar, ZZ r_B_bar, long win_LL){
	long i;
	Mod_p t_B, co_B, temp, temp_1;
	long m = c_B->size();
	vector<Mod_p>* B_temp = new vector<Mod_p>(m+1);
	vector<ZZ>* chal_mult = new vector<ZZ>(m+1);

	//check for correctness of the committed B
	B_temp ->at(0)= c_B0;
	chal_mult->at(0)=1;
	for(i = 0; i<m; i++){
		B_temp->at(i+1)=c_B->at(i);
		chal_mult->at(i+1)=chal->at(i);
	}
	multi_expo::multi_expo_LL(t_B, B_temp, chal_mult, win_LL);

	co_B = Ped.commit_opt(B_bar, r_B_bar);
	//cout<<"B "<<t_B<<" "<<co_B<<endl;
	delete B_temp;
	delete chal_mult;
	if(t_B==co_B){
		return 1;
	}
	return 0;
}

int func_ver::check_a_op(vector<Mod_p>* c_a, vector<Mod_p>* c_Dl, vector<ZZ>* chal, ZZ a_bar, ZZ r_a_bar){
	long i,l;
	Mod_p t_a, co_a, temp, temp_1;
	long m = c_a->size()/2;

	//Check that the random values are used right
	temp = c_a->at(0);
	l=2*m-1;
	for(i = 1; i<=l; i++){
		Mod_p::expo(temp_1,c_a->at(i), chal->at(i-1));
		Mod_p::mult(temp,temp,temp_1);
	}
	t_a = temp;
	co_a = Ped.commit_sw(a_bar, r_a_bar);
	//cout<<"a "<<t_a<<" "<<co_a<<endl;

	if(t_a==co_a & c_a->at(m)==c_Dl->at(m+1)){//both commitments should be com(0,0), so equal
		return 1;
	}
	return 0;
}


int func_ver::check_c_op(vector<vector<Cipher_elg>* >* enc, vector<Cipher_elg>* E, ZZ chal, long omega){
	long i,j;
	ZZ chal_temp;
	ZZ ord = H.get_ord();
	Cipher_elg c, temp;
	vector<ZZ>* v_chal = 0;
	long m = enc->size();
	long n = enc->at(0)->size();

	chal_temp = to_ZZ(1);
	c = Cipher_elg(1,1,H.get_mod());
	v_chal = new vector<ZZ>(n);
 	for(i = 0; i<m;i++){
		for(j = 0; j<n;j++){
			MulMod(chal_temp, chal_temp, chal, ord);
			v_chal->at(j)=chal_temp;
		}
		multi_expo::expo_mult(temp, enc->at(i), v_chal, omega);
		Cipher_elg::mult(c,c,temp);
	}
 	delete v_chal;
 	//cout<<"c "<<E->at(m)<<" "<<c<<endl;
 	if(E->at(m)==c){
 		return 1;
 	}
 	return 0;
}

int func_ver::check_E_op(vector<vector<Cipher_elg>* >* C, vector<Cipher_elg>* E, vector<ZZ>* chal, vector<ZZ>* B_bar, ZZ a_bar, ZZ rho_bar , long omega){
	long i,l;
	Cipher_elg t_E, co_E, temp, temp_1;
	Mod_p gen,t;
	long m = C->size();
	ZZ ord = H.get_ord();
	ZZ mod = H.get_mod();
	long num_b = NumBits(ord);
	l=2*m;
	vector<vector<long>* >* basis_chal = new vector<vector<long>* >(l);

	basis_chal ->at(0)= multi_expo::to_basis(to_ZZ(1), num_b, omega);
	basis_chal->at(1) = multi_expo::to_basis(chal->at(0), num_b, omega);
	for(i=2; i<l; i++){
		basis_chal->at(i) = multi_expo::to_basis(chal->at(i-1), num_b, omega);
	}

	multi_expo::expo_mult(t_E, E, basis_chal, omega);

	vector<ZZ>* chal_temp = new vector<ZZ>(m);
	l=m-1;
	for(i=0; i<l; i++){
		chal_temp->at(i)= chal->at(m-i-2);
	}
	chal_temp->at(l)=1;

	gen = H.get_gen();
	t = Mod_p::expo(gen,a_bar);
	temp = El.encrypt(t, rho_bar);
	multi_expo::expo_mult(temp_1, C, chal_temp, B_bar, omega);

	Cipher_elg::mult(co_E ,temp_1, temp);

	Functions::delete_vector(basis_chal);
	delete(chal_temp);
	//cout<<"E"<<t_E<<" "<<co_E<<endl;
	if(t_E==co_E){
		return 1;
	}

	return 0;

}


void func_ver::fill_vector(vector<ZZ>* t){
	long i,l;
	ZZ temp;
	ZZ ord = H.get_ord();

	l= t->size();
	temp = RandomBnd(ord);
	t->at(0)=temp;
	for(i=1; i<l; i++){
		MulMod(t->at(i),t->at(i-1),temp, ord);
	}
}


void func_ver::fill_x8(vector<ZZ>* chal_x8, vector<vector<long>* >* basis_chal_x8, vector<ZZ>* mul_chal_x8, long omega){
	long i, l;
	ZZ chal;
	ZZ ord = H.get_ord();
	long num_b= NumBits(ord);

	l= chal_x8->size();
	chal = RandomBnd(ord);

	chal_x8->at(0)= chal;
	basis_chal_x8->at(0) = multi_expo::to_basis(to_ZZ(1),num_b, omega);
	basis_chal_x8->at(1) = multi_expo::to_basis(chal_x8->at(0),num_b, omega);

	mul_chal_x8->at(0) =1;
	mul_chal_x8->at(1) =chal_x8->at(0);

	for (i = 1; i<l; i++){
		 MulMod(chal_x8->at(i),chal, chal_x8->at(i-1), ord);
		 basis_chal_x8->at(i+1) = multi_expo::to_basis(chal_x8->at(i),num_b, omega);
		 mul_chal_x8->at(i+1) = chal_x8->at(i);
	}
}

