double rooterfc(double* x, double* par){
    double value = TMath::Erfc((x[0]-par[0])/par[1]);
    //if(value==0) {cout<<x[0]<<endl; exit(0);}
    return value*par[2];
}

void processonefile(string folder, string runnumber){
  gStyle->SetOptFit(1111);

  string head = "/media/sf_H_DRIVE/PROSPECT_LiLS_PSD_data/";
  string rtfilename = head + folder + "/" + runnumber + "_runCs137.root";
  string txtfilename = folder + "_" + runnumber + "_energyScale_varyT.txt";
  TFile* file = new TFile(rtfilename.c_str(),"read");
  fstream fout(txtfilename.c_str(),ios::out);
//  TFile* file1 = new TFile("LiLS_batch2_sample1/run3571741191_shorter.root","read");
  TH1D* hQtot[9];// = (TH1F*)file->Get("hCharge");
  TF1* fErfc[9];
  int T[9]={200,300,400,500,600,700,800,900,1000};
  char name[100];
  for(int i=0;i<9;i++){
    sprintf(name,"hTotCharge_T%d",T[i]);
    hQtot[i] = (TH1D*)file->Get(name);
    hQtot[i]->GetXaxis()->SetRangeUser(0,40);

    sprintf(name,"fErfc%d",T);
    fErfc[i] = new TF1(name,rooterfc,0,50,3);
    fErfc[i]->SetParameter(0, 15.0);
    fErfc[i]->SetParameter(1, 3.0);
    fErfc[i]->SetParameter(2, 160.0);
    hQtot[i]->Fit(name,"SQR","",11,24);
    fErfc[i]->SetLineStyle(2);

    fout<<T[i]<<"\t"<<fErfc[i]->GetParameter(0)<<"\t"<<fErfc[i]->GetParError(0)
              <<"\t"<<fErfc[i]->GetParameter(1)<<"\t"<<fErfc[i]->GetParError(1)
              <<"\t"<<fErfc[i]->GetParameter(2)<<"\t"<<fErfc[i]->GetParError(2)<<endl;

  }
  TCanvas* cQtot = new TCanvas();
  cQtot->Divide(4,2);
  for(int i=0;i<8;i++){
    cQtot->cd(i+1);
    hQtot[i]->Draw();
    hQtot[i]->SetXTitle("Charge (Volts / 50 #Omega*0.2 ns)");
    cQtot->Update();
    TPaveStats* ps=(TPaveStats*)hQtot[i]->GetListOfFunctions()->FindObject("stats");
    ps->SetX1NDC(0.5); ps->SetY1NDC(0.5);
  }

  TCanvas* c1 = new TCanvas();
  hQtot[6]->Draw();

//  TH1F* h1 = (TH1F*)file1->Get("hCharge");
//  h->GetXaxis()->SetRangeUser(0,100);

 // h->Draw();
//  h1->Draw("same"); h1->SetLineColor(kRed);
/*
  TH1F* hdev = new TH1F("hChargeDev","",600,-100,200);
  for(int i=1; i<h->GetNbinsX();i++){
    hdev->SetBinContent(i-1,h->GetBinContent(i)-h->GetBinContent(i-1));
  }

  hdev->Draw();
*/
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
    //if(folder[i]==      "LiLS_batch51_sample1")
      processonefile(folder[i],Cs137Run[i]);
  }
}
