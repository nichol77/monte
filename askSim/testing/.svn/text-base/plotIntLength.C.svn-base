

void plotIntLength() {

   //   TF1 *inty = new TF1("inty",probIntLength,0,1000,1);
   //   inty->SetParameter(0,100);
   //   inty->Draw();
   
   TF1 *proby = new TF1("proby",lengthFromProb,0,1,1);
   proby->SetParameter(0,100);
   proby->Draw();

}


double intLength(double *x, double *par) {
   double L=par[0];
   return TMath::Exp(-1*x[0]/L);
  
}
   

double probIntLength(double *x, double *par) {
   double L=par[0];
   return 1-TMath::Exp(-1*x[0]/L);
  
}


double lengthFromProb(double *p, double *par) {
   double L=par[0];
   return -1*L*TMath::Log(1-p[0]);
}
