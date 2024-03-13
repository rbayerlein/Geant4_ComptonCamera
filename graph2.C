void graph2() {
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,180,700,500);
   c1->SetGrid();
   const Int_t n = 20;
   Double_t x[n], y[n];
   for (Int_t i=0;i<n;i++) {
     x[i] = i*0.1+0.5;
 //    y[i] = 1960*(1-1/(1-(0.511*0.511)/(x[i]+0.511)/(x[i]+0.511))/2.25);
	y[i] = acos (1/(1.5*sqrt(1-0.511/x[i]+0.511)))* 180/3.14159;
     printf(" i %i %f %f \n",i,x[i],y[i]);
   }
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Cherenkov cone opening angle");
   gr->GetXaxis()->SetTitle("Electron Energy [MeV]");
   gr->GetYaxis()->SetTitle(" Cherenkov Angle [degrees]");
   gr->Draw("AC");
   // TCanvas::Update() draws the frame, after which one can change it
   c1->Update();
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
