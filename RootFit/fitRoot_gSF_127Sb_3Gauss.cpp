// Script to fit the 117Sn gamma strength function together with
// 117Sn(gamma,n) data 
// Ann-Cecilie Larsen
// a.c.larsen@fys.uio.no
// 21 Jan 2022
// Last update: 21 Jan 2022


// TOTAL FIT FUNCTION 
// Generalized Lorentzian, one peak - NB!! NOT enhanced GLO! 
// plus pygmy dipole resonance as a Gaussian
// plus low-energy enhancement (upbend)
// 9 parameters althogether

// Function for E1 strength, GLO + pygmy (mainly for plotting, if needed)
double FitFunctionE1(double *x, double *par){
    const double Pi = 3.141592653589793;
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
    const double Pi = 3.141592653589793;
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
    const double Pi = 3.141592653589793;
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
    const double Pi = 3.141592653589793;
    const double factor   = 8.6737E-08;	// const. factor in mb^(-1) MeV^(-2)
    
    // One-component Generalized Lorentzian
    // GDR, as GLO
    // Four parameters:  par[0] = E_r1, par[1] = Gamma_r1, par[2] = sigma_r1, 
    //                   par[3] = temperature of final states 
    double GLo = 0.;
    double GLO_pars[4] = {par[0],par[1],par[2],par[3]};
    GLo = GLO_Function(x,GLO_pars);
    
    // Scissorlike resonance, as a Gaussian
    // Three parameters: par[4] = scaling, par[5] = standard deviation, par[6] = centroid
    double Gauss1 = 0.;
    double PDR_pars1[3] = {par[4],par[5],par[6]};
    Gauss1 = PDR_Function(x,PDR_pars1);

    // Pygmy dipole resonance 1, as a Gaussian
    // Three parameters: par[7] = scaling, par[8] = standard deviation, par[9] = centroid
    double Gauss2 = 0.;
    double PDR_pars2[3] = {par[7],par[8],par[9]};
    Gauss2 = PDR_Function(x,PDR_pars2);
    
    // Pygmy dipole resonance 2, as a Gaussian
    // Three parameters: par[10] = scaling, par[11] = standard deviation, par[12] = centroid
    double Gauss3 = 0.;
    double PDR_pars3[3] = {par[10],par[11],par[12]};
    Gauss3 = PDR_Function(x,PDR_pars3);
     
    // Upbend function, as an exponential
    // Two parameters: par[13] = scaling, par[14] = slope
    double upbend_M1 = 0.;
    double upbend_pars[2] = {par[13],par[14]};
    upbend_M1 = upbend_Function(x,upbend_pars);
   
    double strength_function = GLo + Gauss1 + Gauss2 + Gauss3 + upbend_M1;
    
    return strength_function;
    
}	

void fitRoot_gSF_127Sb_3Gauss(){
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

    // Lepretre data
    //i=0;
    //while(gamma_n_data_file_lepretre){
    //    gamma_n_data_file_lepretre >> x >> y >> z;
    //    energy_lepretre[i] = x; // gamma energy (MeV)
    //    gdr_lepretre[i] = y; // cross-section data already converted to gSF in MeV^(-3)
    //    gdr_lepretre_err[i] = z; // uncertainty
    //    i++;
    //}
    //gamma_n_data_file_lepretre.close();

    // Fultz data
    //i=0;
    //while(gamma_n_data_file_fultz){
    //    gamma_n_data_file_fultz >> x >> y >> z;
    //    energy_fultz[i]  = x; // gamma energy (MeV)
    //    gdr_fultz[i]     = y; // cross-section data already converted to gSF in MeV^(-3)
    //    gdr_fultz_err[i] = z; // uncertainty
    //    i++;
    //}
    //gamma_n_data_file_fultz.close();
    
    // Utsunomiya data
    // This data set is not converted to gSF yet, 
    // will do it while reading in the values from file.
    // Need to skip first lines with text
    //for(i=0;i<13;i++){
    //    getline(gamma_n_data_file_utsunomiya,line);
    //    //cout << line << endl;
    //}
    //i=0;
    //while(gamma_n_data_file_utsunomiya){
    //    gamma_n_data_file_utsunomiya >> x >> y >> z;
    //    energy_utsunomiya[i]  = x;
    //    gdr_utsunomiya[i]     = y*factor/x;
    //    gdr_utsunomiya_err[i] = z*factor/x;
    //    i++;
    //}
    //gamma_n_data_file_utsunomiya.close();


    // Make graphs
    TGraphErrors *strengthexp_oslo = new TGraphErrors(len_oslodata,energy_oslo,strength_oslo,energyerr,strength_oslo_err);
    //TGraphErrors *gdr_lepretre_graph = new TGraphErrors(len_lepretre,energy_lepretre,gdr_lepretre,energyerr,gdr_lepretre_err);
    //TGraphErrors *gdr_fultz_graph = new TGraphErrors(len_fultz,energy_fultz,gdr_fultz,energyerr,gdr_fultz_err);
    //TGraphErrors *gdr_utsunomiya_graph = new TGraphErrors(len_utsunomiya,energy_utsunomiya,gdr_utsunomiya,energyerr,gdr_utsunomiya_err);

    // Make a combined TGraphErrors for the (gamma,n) data, to determine the GDR strength
    //TMultiGraph *GDR_data = new TMultiGraph();
    //GDR_data->Add(gdr_lepretre_graph,"P");
    //GDR_data->Add(gdr_fultz_graph,"P");
    //GDR_data->Add(gdr_utsunomiya_graph,"P");
    
    // Make a combined graph with all data
    TMultiGraph *OCL_and_GDR_data = new TMultiGraph();
    OCL_and_GDR_data->Add(strengthexp_oslo,"P");
    //OCL_and_GDR_data->Add(gdr_lepretre_graph,"P");
    //OCL_and_GDR_data->Add(gdr_fultz_graph,"P");
    //OCL_and_GDR_data->Add(gdr_utsunomiya_graph,"P");
   
   
    // START VALUES for the GDR parameters. 
    // Initializing from 1st fit with just GDR data
    // can give better start values, if it is hard to make the 
    // minimization converge
    double E_r1 = 15.36100; // Centroid in MeV
    double Gamma_r1 = 5.374713;  // Width in MeV
    double sigma_r1 = 283.1516; // Peak cross section in mb
    double temp     = 0.301596; // initial constant temperature of final states (MeV)
    

    // START VALUES Gauss1
    double C_1 = 8.358E-09;
    double sigma_1 = 0.5124;
    double E0_1 = 3.941;

    // START VALUES Gauss2
    double C_2 = 2.8E-07;
    double sigma_2 = 0.87;
    double E0_2 = 6.6;
    
    // START VALUES Gauss3
    double C_3 = 9E-08;
    double sigma_3 = 0.4;
    double E0_3 = 7.8;
    
    // START VALUES upbend
    double constant1_upb = 2.234E-08; // scaling
    double constant2_upb = 0.5; // slope

    
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
    
    //gdr_lepretre_graph->SetMarkerStyle(20);
    //gdr_lepretre_graph->SetMarkerSize(0.8);
    //gdr_lepretre_graph->SetMarkerColor(kAzure+10);
    //gdr_lepretre_graph->SetLineColor(kAzure+10);
    //gdr_lepretre_graph->Draw("P");

    //gdr_fultz_graph->SetMarkerStyle(32);
    //gdr_fultz_graph->SetMarkerSize(1.);
    //gdr_fultz_graph->SetMarkerColor(kMagenta);
    //gdr_fultz_graph->SetLineColor(kMagenta);
    //gdr_fultz_graph->Draw("P");
    
    //gdr_utsunomiya_graph->SetMarkerStyle(33);
    //gdr_utsunomiya_graph->SetMarkerSize(1.4);
    //gdr_utsunomiya_graph->SetMarkerColor(kCyan+2);
    //gdr_utsunomiya_graph->SetLineColor(kCyan+2);
    //gdr_utsunomiya_graph->Draw("P");
    

    // Fit function for the GLO strength (4 parameters)
    // This can first be fit to the (gamma,n) data to get a better estimate for
    // the GDR start parameters
    // Four parameters: par[0] = centroid, par[1] = width, par[2] = peak cross section, par[3] = temperature
    //                              name of new func, func to use, E_start, E_stop, num_parameters
    //TF1 *fit_strength_GLO = new TF1("fit_strength_GLO",GLO_Function,10.,17.,4); 
    //Vector for 4 start parameters (GLO)
    //double par_array_GLO[7] = {E_r1,Gamma_r1,sigma_r1,temp};
    //                       p0     p1       p2    p3     p4      p5      p6  
    //fit_strength_GLO->SetParameters(par_array_GLO);
    // If needed, set limits on the temperature parameter to ease the fit
    //fit_strength_GLO->SetParLimits(3,0.1,1.); 
    // Or fix it if it doesn't behave well
    //fit_strength_GLO->FixParameter(3,0.3); 
    
    //cout << " Fit of GDR part only:" << endl;    
    //GDR_data->Fit(fit_strength_GLO,"RM+");
    
    //fit_strength_GLO->SetLineColor(kBlack);
    //fit_strength_GLO->SetLineWidth(1);
    //fit_strength_GLO->SetLineStyle(1);
    //fit_strength_GLO->Draw("L same"); // Will rather draw the final fit later... 
    
    // Now we have good start parameters for the GLO part.
    // Grab the obtained parameters to be used in the fit 
    // to all data.  
    //E_r1 = fit_strength_GLO->GetParameter(0);
    //Gamma_r1 = fit_strength_GLO->GetParameter(1);
    //sigma_r1 = fit_strength_GLO->GetParameter(2);
    //temp = fit_strength_GLO->GetParameter(3);
    
    // Total fit function with GLO, Gauss1, Gauss2, Gauss3, and upbend. 
    // In total 15 parameters
    TF1 *fit_strength = new TF1("fit_strength",FitFunctionStrength,1.6,16.9,15); 
    
    // Vector for 15 start parameters
    // Pygmy dipole resonance, as a Gaussian
    // Three parameters: par[4] = scaling, par[5] = standard deviation, par[6] = centroid
    double par_array_totalfit[15] = {E_r1,Gamma_r1,sigma_r1,temp,C_1,sigma_1,E0_1,C_2,sigma_2,E0_2,C_3,sigma_3,E0_3,constant1_upb,constant2_upb};
    // In terminal output:           p0     p1       p2      p3  p4   p5      p6  p7    p8     p9  p10   p11   p12     p13            p14            
    fit_strength->SetParameters(par_array_totalfit);
    fit_strength->SetLineColor(kAzure+7);
    // If needed, it's possible to fix parameters from the fit of just the GDR data 
    //fit_strength->FixParameter(0,E_r1); 
    //fit_strength->FixParameter(1,Gamma_r1); 
    //fit_strength->FixParameter(2,sigma_r1); 
    
    // We might also want to set limits on the GDR parameters 
    fit_strength->SetParLimits(0,15.,15.72); 
    fit_strength->SetParLimits(1,4.3497,6.3997);
    fit_strength->SetParLimits(2,254.9266,311.3766);
    fit_strength->SetParLimits(3,0.00,1.3);
    
    //set limits upbend
    fit_strength->SetParLimits(13,1.5E-8,4.7E-8);     
    fit_strength->SetParLimits(14,0.1,1.2); 
    
    cout << endl;      
    cout << " Total fit funtion: " << endl;      
    OCL_and_GDR_data->Fit(fit_strength,"RM+");
        // Also grab the chi^2 and the degrees of freedom to print to the canvas
    TF1 *fitresult = OCL_and_GDR_data->GetFunction("fit_strength");
    double chi2 = fitresult->GetChisquare();
    double degrees_of_freedom = fitresult->GetNDF();
    cout << " The chi^2 value: " << chi2 << endl;
    cout << " Degrees of freedom: " << degrees_of_freedom << endl;
    cout << " Reduced chi^2 = " << chi2/degrees_of_freedom << endl;
    fit_strength->Draw("L same");

 
    // Now we have all the parameters for the fitted gamma strength. 
    // We grab the parameters and make nice functions to draw 
    // together with the data
    E_r1 = fit_strength->GetParameter(0);
    Gamma_r1 = fit_strength->GetParameter(1);
    sigma_r1 = fit_strength->GetParameter(2);
    temp = fit_strength->GetParameter(3);
    C_1 = fit_strength->GetParameter(4);
    sigma_1 = fit_strength->GetParameter(5);
    E0_1 = fit_strength->GetParameter(6);
    C_2 = fit_strength->GetParameter(7);
    sigma_2 = fit_strength->GetParameter(8);
    E0_2 = fit_strength->GetParameter(9);
    C_3 = fit_strength->GetParameter(10);
    sigma_3 = fit_strength->GetParameter(11);
    E0_3 = fit_strength->GetParameter(12);
    constant1_upb = fit_strength->GetParameter(13);
    constant2_upb = fit_strength->GetParameter(14);

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

    // Making a function just for plotting the 1st Gauss
    TF1 *plot_strength_Gauss1 = new TF1("plot_strength_Gauss1",PDR_Function,1.,7.,3); // Gauss (3 parameters)
    plot_strength_Gauss1->FixParameter(0,C_1);
    plot_strength_Gauss1->FixParameter(1,sigma_1);
    plot_strength_Gauss1->FixParameter(2,E0_1);
    plot_strength_Gauss1->SetLineWidth(2);
    plot_strength_Gauss1->SetLineStyle(7);
    plot_strength_Gauss1->SetLineColor(kAzure+2);
    plot_strength_Gauss1->Draw("L same");
    
    // Making a function just for plotting the 2nd Gauss
    TF1 *plot_strength_Gauss2 = new TF1("plot_strength_Gauss2",PDR_Function,2.,16.,3); // Gauss (3 parameters)
    plot_strength_Gauss2->FixParameter(0,C_2);
    plot_strength_Gauss2->FixParameter(1,sigma_2);
    plot_strength_Gauss2->FixParameter(2,E0_2);
    plot_strength_Gauss2->SetLineWidth(2);
    plot_strength_Gauss2->SetLineStyle(7);
    plot_strength_Gauss2->SetLineColor(kAzure+2);
    plot_strength_Gauss2->Draw("L same");

    // Making a function just for plotting the 2nd Gauss
    TF1 *plot_strength_Gauss3 = new TF1("plot_strength_Gauss3",PDR_Function,2.,16.,3); // Gauss (3 parameters)
    plot_strength_Gauss3->FixParameter(0,C_3);
    plot_strength_Gauss3->FixParameter(1,sigma_3);
    plot_strength_Gauss3->FixParameter(2,E0_3);
    plot_strength_Gauss3->SetLineWidth(2);
    plot_strength_Gauss3->SetLineStyle(7);
    plot_strength_Gauss3->SetLineColor(kAzure+2);
    plot_strength_Gauss3->Draw("L same");
    
    // Making a function just for plotting the upbend separately
    TF1 *plot_strength_upbend = new TF1("plot_strength_upbend",upbend_Function,0.1,8.,2); // upbend (2 parameters)
    plot_strength_upbend->FixParameter(0,constant1_upb);
    plot_strength_upbend->FixParameter(1,constant2_upb);
    plot_strength_upbend->SetLineWidth(1);
    plot_strength_upbend->SetLineStyle(2);
    plot_strength_upbend->SetLineColor(kBlack);
    plot_strength_upbend->Draw("L same");

    // Legend for data.         x1   y1   x2  y2   Coordinates in normalized values, x and y from 0 to 1 
    TLegend *leg = new TLegend(0.15,0.70,0.6,0.88);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(strengthexp_oslo," ^{117}Sn(^{3}He,^{3}He'#gamma), OCL ","PE");
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
    leg2->AddEntry(plot_strength_Gauss1," Gauss1 ","L");   
    leg2->AddEntry(plot_strength_Gauss2," Gauss2 ","L"); 
    leg2->AddEntry(plot_strength_Gauss3," Gauss3 ","L"); 
    leg2->AddEntry(plot_strength_upbend," Upbend ","L");   
    leg2->AddEntry(fit_strength," Total fit  ","L");   
    leg2->Draw();


	TLatex t;
	t.SetTextSize(0.05);
	//t.DrawLatex(    19.416,9.080e-07,"^{51}Ti");
    
	c1->Update();
    c1->Print("fitRoot_gSF_117Sn.pdf");

    // Print E1 and M1 strengths to file, to use in the TALYS calculations
    // REMEMBER that the TALYS functions are given in mb/MeV (Goriely's tables)
    // so we must convert with the factor const double factor   = 8.6737E-08;   // const. factor in mb^(-1) MeV^(-2)
    //FILE *E1file, *M1file;  
    
    // We make a function for the total E1 strength, GLO+PDR, 7 parameters
    //TF1 *function_E1 = new TF1("function_E1",FitFunctionE1,0.0,30.,7); 
    // Fix parameters to the ones from the total fit
    //function_E1->FixParameter(0,E_r1); 
    //function_E1->FixParameter(1,Gamma_r1); 
    //function_E1->FixParameter(2,sigma_r1); 
    //function_E1->FixParameter(3,temp); 
    //function_E1->FixParameter(4,scaling1);
    //function_E1->FixParameter(5,stdev1);
    //function_E1->FixParameter(6,E_pyg1);

    //E1file = fopen("E1_gsf_117Sn_TALYSformat.txt","w");
    //fprintf(E1file," Z=  50 A= 117\n");
    //fprintf(E1file,"  U[MeV]  fE1[mb/MeV]\n");
    //double dummy = 0.1;
    //cout << " E1 strength:" << endl;
    //for(i=0;i<300;i++){
    //    fprintf(E1file,"%9.3f%12.3E\n",dummy,(function_E1->Eval(dummy))/factor);
    //    //cout << dummy << " " << (function_E1->Eval(dummy)) << endl;
    //    dummy += 0.1;
    //}
    //fclose(E1file);

    // The M1 part is here believed to be only the upbend. 
    // We could also have included a spin-flip 
    //M1file = fopen("M1_gsf_117Sn_TALYSformat.txt","w");
    //fprintf(M1file," Z=  50 A= 117\n");
    //fprintf(M1file,"  U[MeV]  fM1[mb/MeV]\n");
    //dummy = 0.1;
    //cout << " M1 strength:" << endl;
    //for(i=0;i<300;i++){
    //    fprintf(M1file,"%9.3f%12.3E\n",dummy,(plot_strength_upbend->Eval(dummy))/factor);
    //    //cout << dummy << " " << (plot_strength_upbend->Eval(dummy)) << endl;
    //    dummy += 0.1;
    //}
    //fclose(M1file);
    
    //cout << " Modeled strength functions for 117Sn written in TALYS format. " << endl;

}
