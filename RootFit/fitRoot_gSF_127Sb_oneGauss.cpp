// Script to fit the 127Sb gamma strength function, based on code from Ann-Cecilie Larsen
// Francesco Pogliano

// TOTAL FIT FUNCTION 
// Generalized Lorentzian, one peak - NB!! NOT enhanced GLO! 
// plus pygmy dipole resonance as a Gaussian
// plus low-energy enhancement (upbend)
// 9 parameters althogether

// Function for E1 strength, GLO + pygmy (mainly for plotting, if needed)
double FitFunctionE1(double *x, double *par){
    const double Pi = 3.14159;
    const double factor   = 8.67373E-08;    // const. factor in mb^(-1) MeV^(-2)
    
    // One-component Generalized Lorentzian
    double GLo1 = 0.;
    double Gamma_k1 = 0.;
    double Gamma_k01 = 0.;
    double denominator1 = 0.;
        
    double GLo = 0.;
    double E1strength = 0.;
    
    // Four parameters:  par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = temperature of final states 
    
    // First part
    Gamma_k1 = par[1]*(pow(x[0],2.0) + (2.0*Pi*par[3])*(2.0*Pi*par[3]))/pow(par[0],2.0);
    Gamma_k01 = par[1]*((2.*Pi*par[3])*(2.*Pi*par[3]))/pow(par[0],2.0);	// Gamma_k for E_gamma=0
    denominator1 = ((pow(x[0],2.) - pow(par[0],2.))*(pow(x[0],2.) - pow(par[0],2.))) + pow(x[0],2.)*pow(Gamma_k1,2.);
    GLo1 = factor*par[2]*par[1]*((x[0]*Gamma_k1)/denominator1 + 0.7*Gamma_k01/pow(par[0],3.)); 
    
    GLo = GLo1;

    // Pygmy dipole resonance, as a Gaussian
    // Three parameters: par[4] = scaling, par[5] = standard deviation, par[6] = centroid
    double ExcessStrength = 0.;
    ExcessStrength = par[4]*exp(-pow((x[0]-par[6]),2.)/(2.*pow(par[5],2.)))/(sqrt(2.*Pi)*par[5]);

    E1strength = GLo+ExcessStrength;
    
    return E1strength;
}

// Function for GLO 
double GLO_Function(double *x, double *par){
    const double Pi = 3.14159;
    const double factor   = 8.67373E-08;    // const. factor in mb^(-1) MeV^(-2)
    
    // One-component Generalized Lorentzian
    double GLo1 = 0.;
    double Gamma_k1 = 0.;
    double Gamma_k01 = 0.;
    double denominator1 = 0.;
        
    double GLo = 0.;
    double E1strength = 0.;
    
    // Four parameters:  par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = temperature of final states 
    
    // First part
    Gamma_k1 = par[1]*(pow(x[0],2.0) + (2.0*Pi*par[3])*(2.0*Pi*par[3]))/pow(par[0],2.0);
    Gamma_k01 = par[1]*((2.*Pi*par[3])*(2.*Pi*par[3]))/pow(par[0],2.0); // Gamma_k for E_gamma=0
    denominator1 = ((pow(x[0],2.) - pow(par[0],2.))*(pow(x[0],2.) - pow(par[0],2.))) + pow(x[0],2.)*pow(Gamma_k1,2.);
    GLo1 = factor*par[2]*par[1]*((x[0]*Gamma_k1)/denominator1 + 0.7*Gamma_k01/pow(par[0],3.)); 
    
    GLo = GLo1;

    return GLo;
}
double upbend_Function(double *x, double *par){
     // Upbend function, as an exponential
    // Two parameters: par[0] = scaling, par[1] = slope
    double upbend_M1 = 0.;
    upbend_M1 = par[0]*exp(-par[1]*x[0]);
    
    return upbend_M1;

}

double PDR_Function(double *x, double *par){
    // Pygmy dipole resonance, as a Gaussian
    // Three parameters: par[0] = scaling, par[1] = standard deviation, par[2] = centroid
    const double Pi = 3.14159;
    double ExcessStrength = 0.;
        
    // Gaussian
    ExcessStrength = par[0]*exp(-pow((x[0]-par[2]),2.)/(2.*pow(par[1],2.)))/(sqrt(2.*Pi)*par[1]);

    return ExcessStrength;
}

double SLO_Function(double *x, double *par){// Standard Lorentzian
    const double factor   = 8.67373E-08;    // const. factor in mb^(-1) MeV^(-2)
    double SLO = 0.;
    // 3 parameters: par[0] = E_SLO, par[1] = Gamma_SLO, par[2] = sigma_SLO
    double denominator = (pow(x[0],2.0) - pow(par[0],2.0))*(pow(x[0],2.0) - pow(par[0],2.0)) + pow(par[1],2.0)*pow(x[0],2.0);
    SLO = factor*par[2]*pow(par[1],2.0)*x[0]/denominator;

    return SLO;
}

double FitFunctionStrength(double *x, double *par){ // x = Egamma, par = array of parameters
    const double Pi = 3.14159;
    const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)
    
    // One-component Generalized Lorentzian
    double GLo1 = 0.;
    double Gamma_k1 = 0.;
    double Gamma_k01 = 0.;
    double denominator1 = 0.;
    double GLo = 0.;
    double E1strength = 0.;
    
    // GDR, as GLO
    // Four parameters:  par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = temperature of final states 
    double GLO_pars[4] = {par[0],par[1],par[2],par[3]};
    GLo = GLO_Function(x,GLO_pars);

    // Pygmy dipole resonance, as a Gaussian
    // Three parameters: par[4] = scaling, par[5] = standard deviation, par[6] = centroid
    double ExcessStrength = 0.;
    double PDR_pars[3] = {par[4],par[5],par[6]};
    ExcessStrength = PDR_Function(x,PDR_pars);
     
    // Upbend function, as an exponential
    // Two parameters: par[7] = scaling, par[8] = slope
    double upbend_M1 = 0.;
    double upbend_pars[2] = {par[7],par[8]};
    upbend_M1 = upbend_Function(x,upbend_pars);
    
    // M1 strength
    double strength_M1 = 0.;
    double M1_pars[3] = {par[9],par[10],par[11]};
    strength_M1 = SLO_Function(x,M1_pars);
   
    double strength_function = GLo + ExcessStrength + upbend_M1 + strength_M1;
    
    return strength_function;
    
}	

void fitRoot_gSF_127Sb_oneGauss(){
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetFillColor(1);
	gStyle->SetPadBorderMode(0);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(20000);
	
    // Open files with data
    // Oslo data
    ifstream strengthfile_oslodata("strength.nrm");
    // Photoneutron data converted to gSF, Saclay group 
    //ifstream gamma_n_data_file_lepretre("117Sn_3He3He_117Sn_rec/GDR_Lepretre_Sn117.txt"); 
    // Photoneutron data converted to gSF, Livermore group
    //ifstream gamma_n_data_file_fultz("117Sn_3He3He_117Sn_rec/GDR_Fultz_Sn117.txt"); 
    // Photoneutron data, NewSUBARU group
    //ifstream gamma_n_data_file_utsunomiya("117Sn_3He3He_117Sn_rec/Utsunomiya2009_117Sn_gn.txt");
        
    // Declaration of various stuff: vectors, constants etc.
    const double factor = 8.674E-08; // const. factor in mb^(-1) MeV^(-2)

    const int len_oslodata   = 62; // length + 1 due to while-loop
    //const int len_lepretre   = 46;
    //const int len_fultz      = 106;
    //const int len_utsunomiya = 25;

    double strength_oslo[len_oslodata],strength_oslo_err[len_oslodata],
    energy_oslo[len_oslodata];
    //double energy_lepretre[len_lepretre], gdr_lepretre[len_lepretre], 
    //gdr_lepretre_err[len_lepretre];
    //double energy_fultz[len_fultz], gdr_fultz[len_fultz], gdr_fultz_err[len_fultz];
    //double energy_utsunomiya[len_utsunomiya], gdr_utsunomiya[len_utsunomiya], 
    //gdr_utsunomiya_err[len_utsunomiya];
    // Make large vector with zeros for the energy error (typically not given or ignored)
    double energyerr[200] = {0.0}; 
    
    // Some help variables
    string line;
    line.resize(256);

    int i;
	double x,y,z;
    const double a0 =   -0.7599; // shift (MeV), Oslo data
    const double a1 =   0.1562; // dispersion (MeV/channel), Oslo data


    // Read files and fill vectors.
    // Oslo data
    i=0;
    while(strengthfile_oslodata){
        strengthfile_oslodata >> x;
        if(i<(len_oslodata-1)){
            strength_oslo[i] =   x; // gSF in MeV^(-3)
            energy_oslo[i]   = a0 + (a1*(double)i); // gamma energy in MeV
        }   
        else{strength_oslo_err[i-(len_oslodata-1)] = x;}
        i++;
    }
    strengthfile_oslodata.close();


    // Make graphs
    TGraphErrors *strengthexp_oslo = new TGraphErrors(len_oslodata,energy_oslo,strength_oslo,energyerr,strength_oslo_err);
    
    // Make a combined graph with all data
    TMultiGraph *OCL_and_GDR_data = new TMultiGraph();
    OCL_and_GDR_data->Add(strengthexp_oslo,"P");
   
    // START VALUES for the GDR parameters. 
    // Initializing from 1st fit with just GDR data
    // can give better start values, if it is hard to make the 
    // minimization converge
    double E_r1 = 15.36100; // Centroid in MeV
    double Gamma_r1 = 5.374713;  // Width in MeV
    double sigma_r1 = 283.1516; // Peak cross section in mb
    double temp     = 0.301596; // initial constant temperature of final states (MeV)

    // START VALUES pygmy resonance 1 
    double scaling1 = 3.25E-07; // scaling, starting in the 10^-7 MeV^(-3) ballpark
    double stdev1 = 0.9461; // standard deviation (MeV)
    double E_pyg1 = 6.74; // Pygmy centroid (MeV)

    // START VALUES upbend
    double constant1_upb = 1.863E-08; // scaling
    double constant2_upb = 0.28; // slope
    
    // START VALUES M1 strength
    double E_m1 = 9.3025;
    double gamma_m1 = 3.2984;
    double sigma_m1 = 5.4834;

    
    // Make canvas and histogram for axes, and plot stuff
    TCanvas *c1 = new TCanvas("c1","Gamma-ray strength function",800,700);
    // Histogram for drawing nice axes. 
    // Bins on x axis = 100, start_x = 0.0, stop_x = 18.581 (in MeV)
    // Bins on y axis = 100, start_y = 8.454e-10, stop_x = 7.027e-06  (in MeV^-3)
    TH2F *h = new TH2F("h"," ",100,0.0,   18.581,100,8.454e-10,7.027e-06);

	c1->SetLogy();
	c1->SetLeftMargin(0.14);
    c1->SetBottomMargin(0.13);
	h->GetXaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->CenterTitle();
    h->GetXaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitle("#it{E}_{#gamma} (MeV)");
    h->GetYaxis()->SetTitleFont(42);
    h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->CenterTitle();
	h->GetYaxis()->SetTitleOffset(1.5);
	h->GetYaxis()->SetTitle("#it{f}(#it{E}_{#gamma}) (MeV^{-3})");
	h->Draw();
    
	strengthexp_oslo->SetMarkerStyle(21);
	strengthexp_oslo->SetMarkerSize(0.8);
	strengthexp_oslo->Draw("P");
    
    // Total fit function with GLO, Gaussian pygmy, and upbend. 
    // In total 12 parameters
    TF1 *fit_strength = new TF1("fit_strength",FitFunctionStrength,1.6,16.9,12); 
    
    // Vector for 9 start parameters
    // Pygmy dipole resonance, as a Gaussian
    // Three parameters: par[4] = scaling, par[5] = standard deviation, par[6] = centroid
    double par_array_totalfit[12] = {E_r1,Gamma_r1,sigma_r1,temp,scaling1,stdev1,E_pyg1,constant1_upb,constant2_upb,E_m1,gamma_m1,sigma_m1};
    // In terminal output:           p0     p1       p2     p3     p4      p5     p6         p7        p8         
    fit_strength->SetParameters(par_array_totalfit);
    fit_strength->SetLineColor(kAzure+7);
    // If needed, it's possible to fix parameters from the fit of just the GDR data 
    fit_strength->FixParameter(9,E_m1); 
    fit_strength->FixParameter(10,gamma_m1); 
    fit_strength->FixParameter(11,sigma_m1); 
    
    // We might also want to set limits on the GDR parameters 
    fit_strength->SetParLimits(0,15.,15.72); 
    fit_strength->SetParLimits(1,4.3497,6.3997);
    fit_strength->SetParLimits(2,254.9266,311.3766);
    fit_strength->SetParLimits(3,0.00,1.3);
    
    cout << endl;      
    cout << " Total fit funtion: " << endl;      
    OCL_and_GDR_data->Fit(fit_strength,"RM+");
    fit_strength->Draw("L same");
    TF1 *fitresult = OCL_and_GDR_data->GetFunction("fit_strength");
    double chi2 = fitresult->GetChisquare();
    double degrees_of_freedom = fitresult->GetNDF();
    cout << " The chi^2 value: " << chi2 << endl;
    cout << " Degrees of freedom: " << degrees_of_freedom << endl;
    cout << " Reduced chi^2 = " << chi2/degrees_of_freedom << endl;
    

 
    // Now we have all the parameters for the fitted gamma strength. 
    // We grab the parameters and make nice functions to draw 
    // together with the data
    E_r1 = fit_strength->GetParameter(0);
    Gamma_r1 = fit_strength->GetParameter(1);
    sigma_r1 = fit_strength->GetParameter(2);
    temp = fit_strength->GetParameter(3);
    scaling1 = fit_strength->GetParameter(4);
    stdev1 = fit_strength->GetParameter(5);
    E_pyg1 = fit_strength->GetParameter(6);
    constant1_upb = fit_strength->GetParameter(7);
    constant2_upb = fit_strength->GetParameter(8);

    // Making a function just for plotting the GLO strength for a wider range of energies
    TF1 *plot_strength_GLO = new TF1("plot_strength_GLO",GLO_Function,0.5,25.,4); 
    // Now, fix parameters from the fit above 
    plot_strength_GLO->FixParameter(0,E_r1); 
    plot_strength_GLO->FixParameter(1,Gamma_r1); 
    plot_strength_GLO->FixParameter(2,sigma_r1); 
    plot_strength_GLO->FixParameter(3,temp); 

    plot_strength_GLO->SetLineWidth(1);
    plot_strength_GLO->SetLineColor(kBlack);
    plot_strength_GLO->Draw("L same");

    // Making a function just for plotting the PDR separately
    TF1 *plot_strength_PDR = new TF1("plot_strength_PDR",PDR_Function,2.,16.,3); // PDR (3 parameters)
    plot_strength_PDR->FixParameter(0,scaling1);
    plot_strength_PDR->FixParameter(1,stdev1);
    plot_strength_PDR->FixParameter(2,E_pyg1);
    plot_strength_PDR->SetLineWidth(2);
    plot_strength_PDR->SetLineStyle(7);
    plot_strength_PDR->SetLineColor(kAzure+2);
    plot_strength_PDR->Draw("L same");

    // Making a function just for plotting the upbend separately
    TF1 *plot_strength_upbend = new TF1("plot_strength_upbend",upbend_Function,0.1,8.,2); // upbend (2 parameters)
    plot_strength_upbend->FixParameter(0,constant1_upb);
    plot_strength_upbend->FixParameter(1,constant2_upb);
    plot_strength_upbend->SetLineWidth(1);
    plot_strength_upbend->SetLineStyle(2);
    plot_strength_upbend->SetLineColor(kBlack);
    plot_strength_upbend->Draw("L same");
    
    // Making a function just for plotting the SLO M1 strength for a wider range of energies
    TF1 *plot_strength_SLO = new TF1("plot_strength_SLO",SLO_Function,0.5,25.,3); 
    // Now, fix parameters from the fit above 
    plot_strength_SLO->FixParameter(0,E_m1); 
    plot_strength_SLO->FixParameter(1,gamma_m1); 
    plot_strength_SLO->FixParameter(2,sigma_m1);

    plot_strength_SLO->SetLineWidth(1);
    plot_strength_SLO->SetLineStyle(2);
    plot_strength_SLO->SetLineColor(kBlack);
    plot_strength_SLO->Draw("L same");
   

    // Legend for data.         x1   y1   x2  y2   Coordinates in normalized values, x and y from 0 to 1 
    TLegend *leg = new TLegend(0.15,0.70,0.6,0.88);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(strengthexp_oslo," ^{124}Sn(^{4}He,p'#gamma), OCL ","PE");
    //leg->AddEntry(gdr_lepretre_graph," ^{117}Sn(#gamma,n), Lepretre #it{et al.} ","PE");
    //leg->AddEntry(gdr_fultz_graph," ^{117}Sn(#gamma,n), Fultz #it{et al.} ","PE");  
    //leg->AddEntry(gdr_utsunomiya_graph," ^{117}Sn(#gamma,n), Utsunomiya #it{et al.} ","PE");  
    leg->Draw();

    // Legend for fit functions. x1   y1   x2   y2
    TLegend *leg2 = new TLegend(0.68,0.25,0.85,0.45);
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.03);
    leg2->SetFillColor(0);
    leg2->AddEntry(plot_strength_GLO," GLO  ","L");   
    leg2->AddEntry(plot_strength_PDR," PDR ","L");   
    leg2->AddEntry(plot_strength_upbend," Upbend ","L");   
    leg2->AddEntry(fit_strength," Total fit  ","L");   
    leg2->Draw();
    leg2->AddEntry(plot_strength_SLO," M1 strength  ","L");   

	TLatex t;
	t.SetTextSize(0.05);
	//t.DrawLatex(    19.416,9.080e-07,"^{51}Ti");
    
	c1->Update();
    c1->Print("fitRoot_gSF_127Sb.pdf");

    // Print E1 and M1 strengths to file, to use in the TALYS calculations
    // REMEMBER that the TALYS functions are given in mb/MeV (Goriely's tables)
    // so we must convert with the factor const double factor   = 8.6737E-08;   // const. factor in mb^(-1) MeV^(-2)
    FILE *E1file, *M1file;  
    
    // We make a function for the total E1 strength, GLO+PDR, 7 parameters
    TF1 *function_E1 = new TF1("function_E1",FitFunctionE1,0.0,30.,7); 
    // Fix parameters to the ones from the total fit
    function_E1->FixParameter(0,E_r1); 
    function_E1->FixParameter(1,Gamma_r1); 
    function_E1->FixParameter(2,sigma_r1); 
    function_E1->FixParameter(3,temp); 
    function_E1->FixParameter(4,scaling1);
    function_E1->FixParameter(5,stdev1);
    function_E1->FixParameter(6,E_pyg1);

    E1file = fopen("E1_gsf_117Sn_TALYSformat.txt","w");
    fprintf(E1file," Z=  50 A= 117\n");
    fprintf(E1file,"  U[MeV]  fE1[mb/MeV]\n");
    double dummy = 0.1;
    //cout << " E1 strength:" << endl;
    for(i=0;i<300;i++){
        fprintf(E1file,"%9.3f%12.3E\n",dummy,(function_E1->Eval(dummy))/factor);
        //cout << dummy << " " << (function_E1->Eval(dummy)) << endl;
        dummy += 0.1;
    }
    fclose(E1file);

    // The M1 part is here believed to be only the upbend. 
    // We could also have included a spin-flip 
    M1file = fopen("M1_gsf_127Sb_TALYSformat.txt","w");
    fprintf(M1file," Z=  51 A= 127\n");
    fprintf(M1file,"  U[MeV]  fM1[mb/MeV]\n");
    dummy = 0.1;
    //cout << " M1 strength:" << endl;
    for(i=0;i<300;i++){
        fprintf(M1file,"%9.3f%12.3E\n",dummy,(plot_strength_upbend->Eval(dummy))/factor);
        //cout << dummy << " " << (plot_strength_upbend->Eval(dummy)) << endl;
        dummy += 0.1;
    }
    fclose(M1file);
    
    cout << " Modeled strength functions for 127Sb written in TALYS format. " << endl;

}
