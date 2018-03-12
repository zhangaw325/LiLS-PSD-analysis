double newton(double & y1,double & y2,double & x1,double & x2)
{
    return (y2-y1)/(x2-x1); 
}

void processOneFile(string thefolder, string therun, string theCsrun){
  gStyle->SetOptFit(1111);

  string head = "/media/sf_H_DRIVE/PROSPECT_LiLS_PSD_data/";
  string thename = head + thefolder + "/" + therun + "_runAmBe.root" ;
  string outrootname = thefolder + "_" + therun + "_result_fitLogLikelihood.root";
  string outFOMZname = thefolder + "_" + therun + "_FOM_Z_values.txt";
  string inEnergyscalename = thefolder + "_" + theCsrun + "_energyScale_varyT.txt";

  TFile* rtfile = new TFile(thename.c_str(),"read");
  TFile* rtoutfile = new TFile(outrootname.c_str(),"recreate");
  int T[9]={200, 300, 400, 500, 600, 700, 800, 900, 1000};
  int F[9]={60, 70, 80, 90, 100, 110, 120, 130, 140};
  int chargerange1[7]={5, 15, 25, 35, 45, 55, 65};
  double energyscale[9]; //associated with T array

  fstream fin(inEnergyscalename.c_str(),ios::in);
  fstream fout(outFOMZname.c_str(),ios::out);
  double T0, p0, p1, p2, p3, p4, p5;
  int index = 0;
  while(fin>>T0>>p0>>p1>>p2>>p3>>p4>>p5){
    energyscale[index] = 478.0/p0; //keV per charge unit
    index++;
  }

  TH1D* hPSDvsQtot[9][9][7];

  TGraphErrors* gFOMvsQtot[9][9];
  TGraph* gFOMvsQtot_u[9][9];
  TGraph* gFOMvsQtot_l[9][9];
  TGraphErrors* gZvsQtot[9][9];
  TGraph* gZvsQtot_u[9][9];
  TGraph* gZvsQtot_l[9][9];

  TH2F* hFOMinTFplane = new TH2F("hFOM_TFplane","FOM at 540 keVee",9,60-5,140+5, 9, 150, 1050);
  hFOMinTFplane->SetXTitle("F value");
  hFOMinTFplane->SetYTitle("T value");

  TH2F* hFOMCuminTFplane = new TH2F("hFOMCum_TFplane","Cumulative FOM",9,60-5,140+5, 9, 150, 1050);
  hFOMCuminTFplane->SetXTitle("F value");
  hFOMCuminTFplane->SetYTitle("T value");

  char name[100];
  for(int i=0; i<9; i++){
    //if(i!=6) continue;
    for(int j=0;j<9;j++){
      //if(j!=4) continue;
      sprintf(name,"gFOMvsQtot_T%d_F%d",T[i],F[j]);
      gFOMvsQtot[i][j] = new TGraphErrors();
      gFOMvsQtot[i][j]->SetName(name);
      gFOMvsQtot[i][j]->SetTitle(name);
      sprintf(name,"gFOMvsQtot_u_T%d_F%d",T[i],F[j]);
      gFOMvsQtot_u[i][j] = new TGraph();
      gFOMvsQtot_u[i][j]->SetName(name);
      gFOMvsQtot_u[i][j]->SetTitle(name);
      sprintf(name,"gFOMvsQtot_l_T%d_F%d",T[i],F[j]);
      gFOMvsQtot_l[i][j] = new TGraph();
      gFOMvsQtot_l[i][j]->SetName(name);
      gFOMvsQtot_l[i][j]->SetTitle(name);
      sprintf(name,"gZvsQtot_T%d_F%d",T[i],F[j]);
      gZvsQtot[i][j] = new TGraphErrors();
      gZvsQtot[i][j]->SetName(name);
      gZvsQtot[i][j]->SetTitle(name);
      sprintf(name,"gZvsQtot_u_T%d_F%d",T[i],F[j]);
      gZvsQtot_u[i][j] = new TGraph();
      gZvsQtot_u[i][j]->SetName(name);
      gZvsQtot_u[i][j]->SetTitle(name);
      sprintf(name,"gZvsQtot_l_T%d_F%d",T[i],F[j]);
      gZvsQtot_l[i][j] = new TGraph();
      gZvsQtot_l[i][j]->SetName(name);
      gZvsQtot_l[i][j]->SetTitle(name);
      double fom_sum = 0.0;
      for(int k=0;k<7;k++){
        sprintf(name,"hPSDvsQtot_T%d_F%d_QtotRange_%d_%d",T[i],F[j],chargerange1[k],chargerange1[k]+10);
        hPSDvsQtot[i][j][k] = (TH1D*)rtfile->Get(name);
//        hPSDvsQtot[i][j][k]->Rebin(6);
        double mean = hPSDvsQtot[i][j][k]->GetMean();
        hPSDvsQtot[i][j][k]->Fit("gaus","RQ","",0,mean);
        double c1 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(0);
        double mu1 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(1); 
        double w1 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(2);
        hPSDvsQtot[i][j][k]->Fit("gaus","RQ","",mean,0.4);
        double c2 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(0);
        double mu2 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(1); 
        double w2 = hPSDvsQtot[i][j][k]->GetFunction("gaus")->GetParameter(2);
        TF1* f = new TF1("twogaus","gaus(0)+gaus(3)",0,0.4);
        f->SetParNames("const1","mean1","sigma1","const2","mean2","sigma2");
        f->SetParameters(c1,mu1,w1,c2,mu2,w2);
        //cout<<c1<<"\t"<<mu1<<"\t"<<w1<<"\t"<<c2<<"\t"<<mu2<<"\t"<<w2<<endl;
        TFitResultPtr r = hPSDvsQtot[i][j][k]->Fit("twogaus","SLRQ","",0,0.4);
        TMatrixDSym cov = r->GetCovarianceMatrix();
/*
        r->Write();
        fout<<"Matrix (T "<<T[i]<<"), F "<<F[j]<<")"<<endl;
        for(int iii = 0; iii<6; iii++){
          for(int jjj=0; jjj<6; jjj++){
            fout<<cov[jjj][iii]<<"\t";
          }
          fout<<endl;
        }
*/
        mu1 = f->GetParameter(1); double m1err = f->GetParError(1);
        w1 =  f->GetParameter(2); double w1err = f->GetParError(2);
        mu2 = f->GetParameter(4); double m2err = f->GetParError(4);
        w2 =  f->GetParameter(5); double w2err = f->GetParError(5);
        double Z = mu2-mu1; double Zerr = TMath::Sqrt(m1err**2+m2err**2);
        // 2.3548 = 2*TMath::Sqrt(2.0*TMath::Log(2.0))
        double D = (2.0*TMath::Sqrt(2.0*TMath::Log(2.0)))*TMath::Sqrt((w1)**2+(w2)**2); double Derr = TMath::Sqrt((2.0*(2.3548*w1)*w1err)**2+(2.0*(2.3548*w2)*w2err)**2)/(2.0*D);
        double fom = Z/D;

        double C=2.0*TMath::Sqrt(2.0*TMath::Log(2.0));
        double Q = fom;
        double s2 = w1**2 + w2**2;
        double s = TMath::Sqrt(s2);

        TVector A(6);
        A[0]=0.0; A[1]=-1.0/C/s; A[2]=-1.0*w1*Q/s2;
        A[3]=0.0; A[4]=1.0/C/s;  A[5]=-1.0*w2*Q/s2;
        TVector V0 = A;
        V0 *= cov;
        double myerr = V0*A;
//        double fomerr = fom*TMath::Sqrt((Zerr/Z)**2+(Derr/D)**2);
        double fomerr = TMath::Sqrt(myerr);
        fom_sum += fom;
        gFOMvsQtot[i][j]->SetPoint(k,chargerange1[k]+5,fom);
        gFOMvsQtot[i][j]->SetPointError(k,5,fomerr);
        gFOMvsQtot_u[i][j]->SetPoint(k,chargerange1[k]+5,fom+fomerr);
        gFOMvsQtot_l[i][j]->SetPoint(k,chargerange1[k]+5,fom-fomerr);

        gZvsQtot[i][j]->SetPoint(k,chargerange1[k]+5,Z);
        gZvsQtot[i][j]->SetPointError(k,5,Zerr);
        gZvsQtot_u[i][j]->SetPoint(k,chargerange1[k]+5,Z+Zerr);
        gZvsQtot_l[i][j]->SetPoint(k,chargerange1[k]+5,Z-Zerr);

        hPSDvsQtot[i][j][k]->Write();

        //cout<<fomerr<<"\t"<<TMath::Sqrt(myerr)<<endl;
         
      }
      //prepare to evaluate FOM and error at 540keV 
      double fom_540 = gFOMvsQtot[i][j]->Eval(540.0/energyscale[i]);
      double fom_540e = 0.5*(TMath::Abs(gFOMvsQtot_u[i][j]->Eval(540.0/energyscale[i])-fom_540)+TMath::Abs(gFOMvsQtot_l[i][j]->Eval(540.0/energyscale[i])-fom_540));
      double Z_540 = gZvsQtot[i][j]->Eval(540./energyscale[i]);
      double Z_540e = 0.5*(TMath::Abs(gZvsQtot_u[i][j]->Eval(540.0/energyscale[i])-Z_540)+TMath::Abs(gZvsQtot_l[i][j]->Eval(540.0/energyscale[i])-Z_540));
      fout<<T[i]<<"\t"<<F[j]<<"\t"<<energyscale[i]<<"\t"<<fom_540<<"\t"<<fom_540e<<"\t"<<Z_540<<"\t"<<Z_540e<<endl;

      sprintf(name,"cgFOMvsQtot_T%d_F%d",T[i],F[j]);
      TCanvas* c = new TCanvas(); c->SetName(name);
      gPad->SetGridx(); gPad->SetGridy();
      gFOMvsQtot[i][j]->Draw("AP");
      gFOMvsQtot[i][j]->GetXaxis()->SetTitle("Total charge (volts/50 #Omega*0.2 ns)");
      gFOMvsQtot[i][j]->GetYaxis()->SetTitle("FOM");
      c->Update();
      TPaveText *gtitle = (TPaveText*)gPad->GetPrimitive("title");
      gtitle->SetX1NDC(0.0);
      gPad->Modified();
      gFOMvsQtot[i][j]->GetXaxis()->SetRangeUser(0,80);
      double max_yscale = TMath::MaxElement(gFOMvsQtot[i][j]->GetN(),gFOMvsQtot[i][j]->GetY());
      gFOMvsQtot[i][j]->GetYaxis()->SetRangeUser(0,max_yscale*1.1);
      TGaxis *axis = new TGaxis(0,max_yscale*1.1,80,max_yscale*1.1,0,80.0*energyscale[i],510,"-L");
      axis->SetTitle("Energy scale (keVee)");
      axis->Draw();
      axis->SetLabelFont(gFOMvsQtot[i][j]->GetXaxis()->GetLabelFont());
      axis->SetTitleFont(42); axis->SetTitleSize(0.04); axis->SetLabelSize(0.04);
      c->Write(); c->Close();

      sprintf(name,"cgZvsQtot_T%d_F%d",T[i],F[j]);
      TCanvas* cZ = new TCanvas(); cZ->SetName(name);
      gPad->SetGridx(); gPad->SetGridy();
      gZvsQtot[i][j]->Draw("AP");
      gZvsQtot[i][j]->GetXaxis()->SetTitle("Total charge (volts/50 #Omega*0.2 ns)");
      gZvsQtot[i][j]->GetYaxis()->SetTitle("Z");
      cZ->Update();
      gtitle = (TPaveText*)gPad->GetPrimitive("title");
      gtitle->SetX1NDC(0.0);
      gPad->Modified();
      double min_yscale = TMath::MinElement(gZvsQtot[i][j]->GetN(),gZvsQtot[i][j]->GetY());
      double max_yscale = TMath::MaxElement(gZvsQtot[i][j]->GetN(),gZvsQtot[i][j]->GetY());
      gZvsQtot[i][j]->GetYaxis()->SetRangeUser(min_yscale*0.9,max_yscale*1.1);
      axis = new TGaxis(0,max_yscale*1.1,80,max_yscale*1.1,0,80.0*energyscale[i],510,"-L");
      axis->SetTitle("Energy scale (keVee)");
      axis->Draw();
      axis->SetLabelFont(gFOMvsQtot[i][j]->GetXaxis()->GetLabelFont());
      axis->SetTitleFont(42); axis->SetTitleSize(0.04); axis->SetLabelSize(0.04);
      cZ->Write(); cZ->Close();

      hFOMinTFplane->Fill(F[j], T[i], gFOMvsQtot[i][j]->Eval(540.0/energyscale[i]));
      hFOMCuminTFplane->Fill(F[j], T[i], fom_sum);
    }
  }
  hFOMinTFplane->SetStats(0);
  hFOMinTFplane->Write();
  hFOMCuminTFplane->SetStats(0);
  hFOMCuminTFplane->Write();
  fin.close();
  fout.close();
}

void fitPSD(){
  const int N = 66;
  string folder[N]={
      "P50-1",
      "P50-1",
      "P50-1",
      "P50-1",
      "P50-1",
      "LiLS_batch1_sample1",
      "LiLS_batch1_sample2",
      "LiLS_batch2_sample1",
      "LiLS_batch3_sample1",
      "LiLS_batch4_sample1",
      "LiLS_batch5_sample1",
      "LiLS_batch6_sample1",
      "LiLS_batch7_sample1",
      "LiLS_batch8_sample1",
      "LiLS_batch9_sample1",
      "LiLS_batch10_sample1",
      "LiLS_batch12_sample1",
      "LiLS_batch13_sample1",
      "LiLS_batch14_sample1",
      "LiLS_batch15_sample1",
      "LiLS_batch16_sample1",
      "LiLS_batch17_sample1",
      "LiLS_batch18_sample1",
      "LiLS_batch19_sample1",
      "LiLS_batch20_sample1",
      "LiLS_batch21_sample1",
      "LiLS_batch22_sample1",
      "LiLS_batch23_sample1",
      "LiLS_batch24_sample1",
      "LiLS_batch25_sample1",
      "LiLS_batch25_sample1_repeat",
      "LiLS_batch25_sample2",
      "LiLS_drum12_sample1",
      "LiLS_batch26_sample1",
      "LiLS_batch27_sample1",
      "LiLS_batch28_sample1",
      "LiLS_batch29_sample1",
      "LiLS_batch30_sample1",
      "LiLS_batch31_sample1",
      "LiLS_batch32_sample1",
      "LiLS_batch33_sample1",
      "LiLS_batch34_sample1",
      "LiLS_batch35_sample1",
      "LiLS_batch36_sample1",
      "LiLS_batch37_sample1",
      "LiLS_batch38_sample1",
      "LiLS_batch39_sample1",
      "LiLS_batch40_sample1",
      "LiLS_batch41_sample1",
      "LiLS_batch42_sample1",
      "LiLS_batch43_sample1",
      "LiLS_batch44_sample1",
      "LiLS_batch45_sample1",
      "LiLS_batch47_sample1",
      "LiLS_batch48_sample1",
      "LiLS_batch49_sample1",
      "LiLS_batch50_sample1",
      "LiLS_batch51_sample1",
      "LiLS_batch52_sample1",
      "LiLS_batch53_sample1",
      "LiLS_batch54_sample1",
      "LiLS_batch55_sample1",
      "LiLS_batch56_sample1",
      "LiLS_batch57_sample1",
      "LiLS_batch58_sample1",
      "LiLS_batch59_sample1"
  };
  string AmBeRun[N]={
      "run3568386939",
      "run3568565429",
      "run3568649727",
      "run3569856498",
      "run3570538570",
      "run3571322081",
      "run3571665388",
      "run3571748994",
      "run3571835508",
      "run3571917776",
      "run3571961506",
      "run3572274897",
      "run3572299453",
      "run3572868578",
      "run3572904539",
      "run3573476634",
      "run3574807487",
      "run3574851408",
      "run3574771357",
      "run3574683370",
      "run3575292183",
      "run3575370827",
      "run3575451175",
      "run3575893683",
      "run3575912813",
      "run3576495873",
      "run3577099972",
      "run3577110578",
      "run3577195912",
      "run3577708747",
      "run3577807179",
      "run3579524482",
      "run3582042442",
      "run3577718283",
      "run3577727981",
      "run3577739852",
      "run3579534033",
      "run3579542057",
      "run3579552150",
      "run3579560793",
      "run3579602567",
      "run3579612373",
      "run3580125189",
      "run3580135031",
      "run3580143956",
      "run3581338281",
      "run3581942991",
      "run3581952181",
      "run3581960079",
      "run3581978472",
      "run3582020464",
      "run3582028068",
      "run3582544897",
      "run3583756895",
      "run3583772667",
      "run3583790193",
      "run3583798199",
      "run3583836678",
      "run3583859630",
      "run3584358486",
      "run3584367312",
      "run3584375767",
      "run3586793005",
      "run3586800674",
      "run3594657277",
      "run3594666015"
  };
  string Cs137Run[N]={
      "run3568379453",
      "run3568558558",
      "run3568634069",
      "run3569849534",
      "run3570530846",
      "run3571313172",
      "run3571659302",
      "run3571741191",
      "run3571828361",
      "run3571912967",
      "run3571956769",
      "run3572269892",
      "run3572294635",
      "run3572863130",
      "run3572898975",
      "run3573471805",
      "run3574799917",
      "run3574846553",
      "run3574765657",
      "run3574676766",
      "run3575286093",
      "run3575364948",
      "run3575394729",
      "run3575889206",
      "run3575902045",
      "run3576493547",
      "run3577095349",
      "run3577108138",
      "run3577133217",
      "run3577704489",
      "run3577804902",
      "run3579520791",
      "run3582040946",
      "run3577716724",
      "run3577725850",
      "run3577737823",
      "run3579532323",
      "run3579540579",
      "run3579550074",
      "run3579558576",
      "run3579600617",
      "run3579608948",
      "run3580123381",
      "run3580132658",
      "run3580141459",
      "run3581333443",
      "run3581939990",
      "run3581950320",
      "run3581958497",
      "run3581975415",
      "run3582018679",
      "run3582026572",
      "run3582541957",
      "run3583754608",
      "run3583770150",
      "run3583788243",
      "run3583796474",
      "run3583834257",
      "run3583855739",
      "run3584356009",
      "run3584365793",
      "run3584373409",
      "run3586791220",
      "run3586799184",
      "run3594654725",
      "run3594663218"
  };

  for(int i=0; i<N; i++){
    if(folder[i]==      "LiLS_batch51_sample1")
      processOneFile(folder[i],AmBeRun[i],Cs137Run[i]);
  }
}
