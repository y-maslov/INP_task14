#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TVirtualFFT.h>
#include <vector>


Int_t bins = 500;
Double_t gaus_linear(Double_t *x, Double_t *a)
{
    return a[0] * TMath::Gaus(x[0], a[1], a[2], true) + a[3] + a[4] * x[0];
}


void task14()
{       

    Double_t x[bins + 1]; 
    auto c = new TCanvas("Canvas", "Canvas", 1920, 1020);
    c->Divide(5, 1);
    

    auto h1 = new TH1F("hist", "hist", bins, 0, 500); 
    


    std::ifstream file("dataFFT.dat");
    Float_t x_it = 0;
    Int_t cnt = 0;
    while(file >> x_it)
    {
        x[cnt] = x_it;
        cnt ++;
    }
    file.close();
    for (Int_t i = 0; i < bins; i ++) h1->SetBinContent(i + 1, x[i]);
   
    
    c->cd(1);
    h1->GetXaxis()->SetTitle("ns");
    h1->Draw();
    

    c->cd(2);
    
    TVirtualFFT::SetTransform(0); 
    TH1 *hm = 0;
    hm = h1->FFT(hm, "MAG");
    hm->SetName("spectrum");
    hm->SetTitle("spectrum");
    
  

    hm->GetXaxis()->SetRange(0, (hm->GetEntries()) / 2);
    c->cd(2)->SetLogy();
    hm->Draw();
     
    auto hm_filtered = (TH1F*)hm->Clone();
    hm_filtered->SetName("spectrum filtered");
    hm_filtered->SetTitle("spectrum filtered");
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    Double_t *re_full = new Double_t[bins];
    Double_t *im_full = new Double_t[bins];
    fft->GetPointsComplex(re_full,im_full);
   

    for (Int_t i = 0; i < bins; i ++) 
    {
        if (hm_filtered->GetBinContent(i + 1) < 20) 
        {
            hm_filtered->SetBinContent(i + 1, 0);
            im_full[i] = 0;
            re_full[i] = 0;
        }
       
    }
    hm_filtered->GetXaxis()->SetRange(0, (hm_filtered->GetEntries()) / 2);
    c->cd(3)->SetLogy();
    hm_filtered->Draw();
    TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &bins, "C2R");
    fft_back->SetPointsComplex(re_full,im_full);
    fft_back->Transform();


    c->cd(4);
    TH1 *h2 = 0;
    h2 = TH1::TransformHisto(fft_back,h2,"MAG");
    h2->SetName("signal filtered");
    h2->SetTitle("signal filtered");
  
    h2->Scale(1. / bins);

    TF1 *fit_func = new TF1("fit", *gaus_linear, 0, bins, 6);
    fit_func->SetParameter(0,   1.5);
    fit_func->SetParameter(1,  250);
    fit_func->SetParameter(2, 50);
    fit_func->SetParameter(3, 0);
    fit_func->SetParameter(3, 0);
  
    auto fit = h2->Fit("fit");
    h2->Draw();

    c->cd(5);
    Double_t ampl = fit_func->GetParameter(0);
    Double_t mean = fit_func->GetParameter(1);
    Double_t sigma = fit_func->GetParameter(2);
    TF1 *gaus = new TF1("gaus", "[0] * TMath::Gaus(x, [1], [2], true)", 0, bins);
    gaus->SetParameter(0, ampl);
    gaus->SetParameter(1, mean);
    gaus->SetParameter(2, sigma);
    gaus->SetName("gaus");
    gaus->SetTitle("gaus");
    gaus->GetXaxis()->SetTitle("ns");
    gaus->Draw();

    std::cout << "0.2 ampl time: " << gaus->GetX(0.2 * gaus->GetMaximum()) <<std::endl;
}


