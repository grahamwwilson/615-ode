/// Get in memory an histogram from a root file and fit a user defined function.
/// Note that a user defined function must always be defined
/// as in this example:
///  - first parameter: array of variables (in this example only 1-dimension)
///  - second parameter: array of parameters
/// Note also that in case of user defined functions, one must set
/// an initial value for each parameter.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Rene Brun

Double_t fitf(Double_t *x, Double_t *par)
{
   Double_t t=x[0];     // time in s
   Double_t arg = 0.0;
   Double_t PIE = 4.0*atan(1.0);
   Double_t Omega = 2.0*PIE/par[2]; // par[2] is the period
   arg = Omega*t - par[3];
//
//   if (par[2] != 0) arg = (x[0] - par[1])/par[2];

   Double_t fitval = par[0] + par[1]*TMath::Cos(arg);
   return fitval;
}
void myfit()
{

   TFile *f = new TFile("rk4_0p01.root");
   if (!f) return;

   TCanvas *c1 = new TCanvas("c1","the fit canvas",500,500);

   TH1D *hx = (TH1D*)f->Get("homega");

// Creates a Root function based on function fitf above
   TF1 *func = new TF1("fitf",fitf,0.0,10.0,4);

// Sets initial values and parameter names
   func->SetParameters(0.0,10.0,2.8,0.0);
   func->SetParNames("Mean","Amplitude","Period","Phase");

// Prior to fitting set the bin-by-bin errors appropriately for the histogram
   for (Int_t binx=1; binx<=hx->GetNbinsX(); binx++){
      Int_t global_bin = hx->GetBin(binx);
      hx->SetBinError(global_bin,0.001);
   }

   func->FixParameter(0,0.0);
//   func->FixParameter(2,2.8);
//   func->FixParameter(3,0.0);

// Fit histogram in range defined by function
   hx->Fit(func,"r");

// Gets integral of function between fit limits
//   printf("Integral of function = %g\n",func->Integral(-2,2));

   hx->Draw("hist");
   func->Draw("same");

}
