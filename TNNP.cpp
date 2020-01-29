/*
 * This header implements the ten Tusscher 2006 cardiac cell.
 * K. H. W. J. ten Tusscher and A. V. Panfilov, “Alternans 
 * and spiral breakup in a human ventricular tissue model,” 
 * Am. J. Physiol. Heart Circ. Physiol.
 *
 * Author (Original Version): KHJW ten Tusscher
 * Author (This C++ version): Dominic Whittaker <dominic.whittaker@nottingham.ac.uk>
 * Date Modified			: 29th January, 2020
 *
 * Compile with: g++ TNNP.cpp -o ttcell (on Linux)
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#define OUTFREQ 50

//===========
// Cell types
//===========
enum TypeCell { EPI, MCELL, ENDO };


//==================
// ten Tusscher Cell
//==================
class TTCell {
public:	
	//-------------
	// Constructors
	//-------------
	TTCell(TypeCell celltype = EPI, double epiMidRatio = 1.5, double HT = 0.02) :
	m_celltype(celltype), m_epiMidRatio(epiMidRatio), m_HT(HT)
	{		
		Ko=5.4; 
		Cao=2.0; 
		Nao=140.0; 
		Vc=0.016404; 
		Vsr=0.001094; 
		Vss=0.00005468; 
		Bufc=0.2; 
		Kbufc=0.001; 
		Bufsr=10.; 
		Kbufsr=0.3; 
		Bufss=0.4; 
		Kbufss=0.00025; 
		Vmaxup=0.006375; 
		Kup=0.00025; 
		Vrel=0.102;//*/40.8; 
		k1_=0.15; 
		k2_=0.045; 
		k3=0.060; 
		k4=0.005;//*/0.000015; 
		EC=1.5; 
		maxsr=2.5; 
		minsr=1.; 
		Vleak=0.00036; 
		Vxfer=0.0038; 
		R=8314.472; 
		F=96485.3415; 
		T=310.0; 
		RTONF=(R*T)/F;
		CAPACITANCE=0.185; 
		pKNa=0.03; 
		
		// SHORT AP Gkr*5 Gks*5
		// BASE AP Gks*0.5 Gto*0.5
		if (m_celltype == EPI) {
			Gks=1*0.392; 
			Gkr=1*0.2448;  
		}
		else if (m_celltype == ENDO) {
			Gks=1*0.392;
			Gkr=1*0.153;  
		}
		else if (m_celltype == MCELL) {
			Gks=1*0.098;
			Gkr=1*0.153;   
		}
		
		GK1=5.405; 
		//Parameters for Ito
		if (m_celltype == EPI)
			Gto=1*0.294; 
		else if (m_celltype == ENDO)
			Gto=1*0.073;
		else if (m_celltype == MCELL)
			Gto=1*0.294;
		
		GNa=14.838; 
		GbNa=0.00029; 
		KmK=1.0; 
		KmNa=40.0; 
		knak=2.724;  
		GCaL=0.00003980;  
		GbCa=0.000592; 
		knaca=1000;  
		KmNai=87.5; 
		KmCa=1.38; 
		ksat=0.1; 
		n=0.35; 
		GpCa=0.1238; 
		KpCa=0.0005; 
		GpK=0.0146; 
		
		svolt=-86.2; 
		Cai=0.00007; 
		CaSR=1.3;
		CaSS=0.00007;
		Nai=7.67;
		Ki=138.3;
		
		stimduration=1.;
		stimstrength=-52;
		tbegin=0;
		tend=tbegin+stimduration;
		counter=1;
		count = 0;
		dia=5000;
		basicperiod=1000.0;
		basicapd=274;
		repeats=10;
		
		inverseVcF2=1/(2*Vc*F);
		inverseVcF=1./(Vc*F);
		inversevssF2=1/(2*Vss*F);
		
		sm = 0.0;
		sh = 0.75;
		sj = 0.75;
		sxr1 = 0.0;
		sxr2 = 1.0;
		sxs = 0.0; 
		sr = 0.0;
		ss = 1.0;  
		sd = 0.0;
		sf = 1.0;
		sf2 = 1.0;
		sfcass = 1.0;
		sRR = 1.0;
		sOO = 0.0;
		svolt = -86.2;
		Cai;
		CaSR;
		CaSS;
		Nai;
		Ki;
		sItot = 0.0;

		// APD Parameters
		APD_flag = 0;
		Vo       = svolt;    
		vdot     = 0;
		vdot_old = 0;
		vdot_max = 0;
		t_vdot_max = 0;

	}
	
	//-----------
	// Destructor
	//-----------
	~TTCell() {}		
	
	void CurrentDependencies()
	{
		Ek=RTONF*(log((Ko/Ki)));
		Ena=RTONF*(log((Nao/Nai)));
		Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
		Eca=0.5*RTONF*(log((Cao/Cai)));
		Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
		Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
			 exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
		rec_iK1=Ak1/(Ak1+Bk1);
		rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
		rec_ipK=1./(1.+exp((25-svolt)/5.98));
	}
	
	void CalculateINa() 
	{
		INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
	}
	
	void CalculateICaL()
	{
		ICaL=GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
		(0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
	}
	
	void CalculateIto()
	{
		Ito=Gto*sr*ss*(svolt-Ek);
	}
	
	void CalculateTTIKr()
	{
		IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
	}
	
	void CalculateIKs()
	{
		IKs=Gks*sxs*sxs*(svolt-Eks);
	}
	
	void CalculateIK1()
	{
		IK1=GK1*rec_iK1*(svolt-Ek);
	}
	
	void CalculateINaCa() {
		INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
		(1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
		(exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
		 exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
	}
	
	void CalculateINaK()
	{
		INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
	}
	
	void CalculateIpCa()
	{
		IpCa=GpCa*Cai/(KpCa+Cai);
	}
	
	void CalculateIpK()
	{
		IpK=GpK*rec_ipK*(svolt-Ek);
	}
	
	void CalculateIbNa() 
	{
		IbNa=GbNa*(svolt-Ena);
	}
	
	void CalculateIbCa()
	{
		IbCa=GbCa*(svolt-Eca);
	}

	void CalculateItot()
	{
		sItot = IKr + IKs + IK1 + Ito +	INa + IbNa + ICaL +	IbCa + INaK + INaCa + IpCa + IpK + Istim;
	}
	
	void UpdateConcentrations()
	{
		kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR))); 
		k1=k1_/kCaSR;
		k2=k2_*kCaSR;
		dRR=k4*(1-sRR)-k2*CaSS*sRR;
		sRR+=m_HT*dRR;
		sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
		Irel=Vrel*sOO*(CaSR-CaSS);
		
		Ileak=Vleak*(CaSR-Cai);
		Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
		Ixfer=Vxfer*(CaSS-Cai);
		
		
		CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
		dCaSR=m_HT*(Iup-Irel-Ileak);
		bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
		cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
		CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
		
		
		CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
		dCaSS=m_HT*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
		bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
		ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
		CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
		
		
		CaBuf=Bufc*Cai/(Cai+Kbufc);
		dCai=m_HT*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
		bc=Bufc-CaBuf-dCai-Cai+Kbufc;
		cc=Kbufc*(CaBuf+dCai+Cai);
		Cai=(sqrt(bc*bc+4*cc)-bc)/2;
        
		
		dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
		Nai+=m_HT*dNai;
		
		dKi=-(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
		Ki+=m_HT*dKi;
	}
		
	void ComputeSteadyStateValues()
	{
		AM=1./(1.+exp((-60.-svolt)/5.));
		BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
		TAU_M=AM*BM;
		M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
		if (svolt>=-40.)
		{
			AH_1=0.; 
			BH_1=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
			TAU_H= 1.0/(AH_1+BH_1);
		}
		else
		{
			AH_2=(0.057*exp(-(svolt+80.)/6.8));
			BH_2=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
			TAU_H=1.0/(AH_2+BH_2);
		}
		H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
		if(svolt>=-40.)
		{
			AJ_1=0.;      
			BJ_1=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
			TAU_J= 1.0/(AJ_1+BJ_1);
		}
		else
		{
			AJ_2=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
				   exp(-0.04391*svolt))*(svolt+37.78)/
				  (1.+exp(0.311*(svolt+79.23))));    
			BJ_2=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
			TAU_J= 1.0/(AJ_2+BJ_2);
		}
		J_INF=H_INF;
		
		Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
		axr1=450./(1.+exp((-45.-svolt)/10.));
		bxr1=6./(1.+exp((svolt-(-30.))/11.5));
		TAU_Xr1=axr1*bxr1;
		Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
		axr2=3./(1.+exp((-60.-svolt)/20.));
		bxr2=1.12/(1.+exp((svolt-60.)/20.));
		TAU_Xr2=axr2*bxr2;
		
		Xs_INF=1./(1.+exp((-5.-svolt)/14.));
		Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
		Bxs=(1./(1.+exp((svolt-35.)/15.)));
		TAU_Xs=Axs*Bxs+80;
		
		if (m_celltype == EPI) {
			R_INF=1./(1.+exp((20-svolt)/6.));
			S_INF=1./(1.+exp((svolt+20)/5.));
			TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
			TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
		}
		else if (m_celltype == ENDO) {
			R_INF=1./(1.+exp((20-svolt)/6.));
			S_INF=1./(1.+exp((svolt+28)/5.));
			TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
			TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
		}
		else {
			R_INF=1./(1.+exp((20-svolt)/6.));
			S_INF=1./(1.+exp((svolt+20)/5.));
			TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
			TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
		}		
		
		D_INF=1./(1.+exp((-8-svolt)/7.5));
		Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
		Bd=1.4/(1.+exp((svolt+5)/5));
		Cd=1./(1.+exp((50-svolt)/20));
		TAU_D=Ad*Bd+Cd;
		F_INF=1./(1.+exp((svolt+20)/7));
		Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
		Bf=200./(1+exp((13-svolt)/10.));
		Cf=(180./(1+exp((svolt+30)/10)))+20;
		TAU_F=Af+Bf+Cf;
		F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
		Af2=600*exp(-(svolt+25)*(svolt+25)/170);
		Bf2=31/(1.+exp((25-svolt)/10));
		Cf2=16/(1.+exp((svolt+30)/10));
		TAU_F2=Af2+Bf2+Cf2;
		FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
		TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;

	}
	
	void UpdateGates()
	{
		sm = M_INF-(M_INF-sm)*exp(-m_HT/TAU_M);
		sh = H_INF-(H_INF-sh)*exp(-m_HT/TAU_H);
		sj = J_INF-(J_INF-sj)*exp(-m_HT/TAU_J);
		sxr1 = Xr1_INF-(Xr1_INF-sxr1)*exp(-m_HT/TAU_Xr1);
		sxr2 = Xr2_INF-(Xr2_INF-sxr2)*exp(-m_HT/TAU_Xr2);
		sxs = Xs_INF-(Xs_INF-sxs)*exp(-m_HT/TAU_Xs);
		ss= S_INF-(S_INF-ss)*exp(-m_HT/TAU_S);
		sr= R_INF-(R_INF-sr)*exp(-m_HT/TAU_R);
		sd = D_INF-(D_INF-sd)*exp(-m_HT/TAU_D); 
		sf =F_INF-(F_INF-sf)*exp(-m_HT/TAU_F); 
		sf2 =F2_INF-(F2_INF-sf2)*exp(-m_HT/TAU_F2); 
		sfcass =FCaSS_INF-(FCaSS_INF-sfcass)*exp(-m_HT/TAU_FCaSS);
	}
	
	void UpdateVoltage()
	{
		svolt= svolt + m_HT*(-sItot);
	}
	
	void dVdt_APD( double time )
    {
        double APA = Vmax-Vrest;
        vdot_old = vdot;
        vdot	 = (svolt-svolto)/m_HT;
        if ( (!apd_flag) and (svolt > -40) and (vdot < vdot_old) ) {
            vdot_max	= vdot_old;
            t_vdot_max	= time - m_HT;
            apd_flag	= true;
        }
        if ( apd_flag and (svolt < 0.9*(Vmax - APA))) {
            APD		 = time - t_vdot_max;
            apd_flag = false;
        }
	}

    void setBeats( std::size_t beats ) { repeats = beats; }
    
    void setBCL( double bcl ) { basicperiod = bcl; }
    
	//=========
	// Run Cell
	//=========
	void operator()( std::ofstream &output_file, std::ofstream &IC_file )
    {
        double t = 0;   // Time
        dVdtmax = 0;
		std::size_t total_simulation_steps = static_cast<std::size_t>( repeats * basicperiod/m_HT );

		for ( std::size_t step = 0; step < total_simulation_steps; ++step )
        {
			
			// if (svolt < Vrest) Vrest = svolt;

            if (t >= tbegin and t <= tend)
				Istim = stimstrength;
            if (t > tend) {
                Istim = 0.0;
                tbegin = tbegin + basicperiod;
                tend = tbegin + stimduration;
                ++counter;
                
                std::cout << "Stimulus " << counter << " applied, Time = " << tbegin << std::endl;
            }
			
			//===============
			CurrentDependencies();		
			CalculateINa();
			CalculateICaL();
			CalculateIto();
			CalculateTTIKr();		
			CalculateIKs();
			CalculateIK1();
			CalculateINaCa();
			CalculateINaK();
			CalculateIpCa();
			CalculateIpK();
			CalculateIbNa();
			CalculateIbCa();
			CalculateItot();
			
			UpdateConcentrations();

			// Write Currents To File
            if ( step >= ( static_cast<std::size_t>((repeats-1)*basicperiod/m_HT) - 50/m_HT ) and count % OUTFREQ == 0 ) { // Output last beat
				output_file << t        << " "//std::setw(15)
                            << svolt    << " "//std::setw(15)
                            << IKr      << " "//std::setw(15)
                            << INa      << " "//std::setw(15)
                            << ICaL     << " "//std::setw(15)
                            << Ito      << " "//std::setw(15)
                            << std::endl;
            }
				
			svolto = svolt;
			ComputeSteadyStateValues();
			UpdateGates();
			UpdateVoltage();			
			//===============

			if (step > total_simulation_steps/2) {
				if (svolt < Vrest) Vrest = svolt;
				if (svolt > Vmax) Vmax = svolt;
				dVdt = (svolt - svolto)/m_HT;
				if (dVdt > dVdtmax) dVdtmax = dVdt;
				// APD Calculation
				dVdt_APD( t );
			}
			Vo = svolt;
			t += m_HT;
			count ++;		

		} // END: step loop
		std::cout << "MUV: " << dVdtmax << " V/s" << std::endl;
		std::cout << "APD: " << APD << " ms" << std::endl;

		IC_file << Cai 		<< "\n"
				<< CaSR 	<< "\n"
				<< CaSS 	<< "\n"
				<< sd 		<< "\n"
				<< sf 		<< "\n"
				<< sf2 		<< "\n"
				<< sfcass 	<< "\n"
				<< sh 		<< "\n"
				<< sj 		<< "\n"
				<< Ki 		<< "\n"
				<< sm 		<< "\n"
				<< Nai 		<< "\n"
				<< sr 		<< "\n"
				<< ss 		<< "\n"
				<< svolt 	<< "\n"
				<< sxr1 	<< "\n"
				<< sxr2 	<< "\n"
				<< sxs 		<< "\n"
				<< sRR 		<< "\n"
				<< sOO 		<< "\n"
				<< std::endl;

	}
	
	// private:
	double m_epiMidRatio;
	TypeCell m_celltype;
	
	//External concentrations
	double Ko;
	double Cao;
	double Nao;
	
	//Intracellular volumes
	double Vc;
	double Vsr;
	double Vss;
	
	//Calcium buffering dynamics
	double Bufc;
	double Kbufc;
	double Bufsr;
	double Kbufsr;
	double Bufss;
	double Kbufss;
	
	//Intracellular calcium flux dynamics
	double Vmaxup;
	double Kup;
	double Vrel;
	double k1_;
	double k2_;
	double k3;
	double k4;
	double EC;
	double maxsr;
	double minsr;
	double Vleak;
	double Vxfer;	
	
	//Constants
	double R;
	double F;
	double T;
	double RTONF;
	
	//Cellular capacitance         
	double CAPACITANCE;
	
	//Parameters for currents
	//Parameters for IKr
	double Gkr;
	
	//Parameters for Iks
	double pKNa;
	double Gks;
	
	//Parameters for Ik1
	double GK1;
	
	//Parameters for Ito
	double Gto;
	
	//Parameters for INa
	double GNa;
	
	//Parameters for IbNa
	double GbNa;
	
	//Parameters for INaK
	double KmK;
	double KmNa;
	double knak;
	
	//Parameters for ICaL
	double GCaL;
	
	//Parameters for IbCa
	double GbCa;
	
	//Parameters for INaCa
	double knaca;
	double KmNai;
	double KmCa;
	double ksat;
	double n;
	
	//Parameters for IpCa
	double GpCa;
	double KpCa;
	
	//Parameters for IpK;
	double GpK;
	
	
	//==========================
	// PARAMETER FOR INTEGRATION
	//==========================
	double  m_HT; // time step
	
	//==================================
	// PARAMETERS FOR INITIAL CONDITIONS 
	//==================================
	//Initial values of state variables
	double svolt, svolto;
	double dVdtmax, dVdt;
	double Cai;
	double CaSR;
	double CaSS;
	double Nai;
	double Ki;
	
	//=====================================
	// PARAMETERS FOR STIMULATION PROTOCOLS 
	//=====================================	
	double stimduration;
	double stimstrength;
	double tbegin;
	double tend;
	int counter, count;
	double dia;
	double basicperiod;
	double basicapd;
	int repeats;
	double Istim;
	
	double IKr;
	double IKs;
	double IK1;
	double Ito;
	double INa;
	double IbNa;
	double ICaL;
	double IbCa;
	double INaCa;
	double IpCa;
	double IpK;
	double INaK;
	double Irel;
	double Ileak;
	double Iup;
	double Ixfer;
	double k1;
	double k2;
	double kCaSR;
	
	double dNai;
	double dKi;
	double dCai;
	double dCaSR;
	double dCaSS;
	double dRR;
	
	double Ek;
	double Ena;
	double Eks;
	double Eca;
	double CaCSQN;
	double bjsr;
	double cjsr;
	double CaSSBuf;
	double bcss;
	double ccss;
	double CaBuf;
	double bc;
	double cc;
	double Ak1;
	double Bk1;
	double rec_iK1;
	double rec_ipK;
	double rec_iNaK;
	double AM;
	double BM;
	double AH_1;
	double BH_1;
	double AH_2;
	double BH_2;
	double AJ_1;
	double BJ_1;
	double AJ_2;
	double BJ_2;
	double M_INF;
	double H_INF;
	double J_INF;
	double TAU_M;
	double TAU_H;
	double TAU_J;
	double axr1;
	double bxr1;
	double axr2;
	double bxr2;
	double Xr1_INF;
	double Xr2_INF;
	double TAU_Xr1;
	double TAU_Xr2;
	double Axs;
	double Bxs;
	double Xs_INF;
	double TAU_Xs;
	double R_INF;
	double TAU_R;
	double S_INF;
	double TAU_S;
	double Ad;
	double Bd;
	double Cd;
	double Af;
	double Bf;
	double Cf;
	double Af2;
	double Bf2;
	double Cf2;
	double TAU_D;
	double D_INF;
	double TAU_F;
	double F_INF;
	double TAU_F2;
	double F2_INF;
	double TAU_FCaSS;
	double FCaSS_INF;
	
	double inverseVcF2;
	double inverseVcF;
	double inversevssF2;
	
	double sm;
	double sh;
	double sj;
	double sxr1;
	double sxr2;
	double sxs; 
	double ss;  
	double sr;
	double sd;
	double sf;
	double sf2;
	double sfcass;
	double sRR;
	double sOO;
	double sItot;

	// APD Parameters
    bool apd_flag;
	double vdot, vdot_old, Vo, APD_flag, vdot_max, t_vdot_max, Vrest, Vmax, APD;
};



int main()
{
	std::ofstream output_file( "APs/ENDO.txt" );                        // output filename

    TypeCell celltype       = ENDO;                                 // Type of cell { EPI, MCELL, ENDO }
	double epi_mid_ratio    = 1.5;                                  // EPi:MID ratio
    double time_step        = 0.02;                                 // Time step
    
    TTCell cell( celltype, epi_mid_ratio, time_step );              // Create cell
    
    cell.setBeats( 100 );                                           // Set number of beats
    cell.setBCL( 1000.0 );                                          // Set basic cycle length
    std::ofstream IC_file( "ICs.txt");								// Initial conditions for tissue simulations

    //=========
    // Run cell
    //=========
	cell( output_file, IC_file );
    
    // int plot = std::system( "python plot.py" );                     // Plot results
}