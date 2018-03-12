void draw_psd_histograms(){
    //fstream fin("PSD_data.txt",ios::in);
    string runname;
    float energycal;
    float fom, fomerr;
    float Z, Zerr;
    float fast, total;
    string date1, time1, date2, time2; 
    string filename1, filename2;
    string jobtime;

    const int n = 66;
    const char *people[n] = {"P50-m1","P50-m2","P50-m3","P50-m4","P50-m5", 
                             "B1-S1","B1-S2","B2-S1","B3-S1","B4-S1",
                             "B5-S1","B6-S1","B7-S1","B8-S1","B9-S1",
                             "B10-S1","B12-S1","B13-S1","B14-S1","B15-S1",
                             "B16-S1","B17-S1","B18-S1","B19-S1","B20-S1",
                             "B21-S1","B22-S1","B23-S1","B24-S1","B25-S1",
                             "B25-S1r","B25-S2","D12-S1","B26-S1","B27-S1",
                             "B28-S1","B29-S1","B30-S1","B31-S1","B32-S1",
                             "B33-S1","B34-S1","B35-S1","B36-S1","B37-S1",
                             "B38-S1","B39-S1","B40-S1","B41-S1","B42-S1",
                             "B43-S1","B44-S1","B45-S1","B47-S1","B48-S1",
                             "B49-S1","B50-S1","B51-S1","B52-S1","B53-S1",
                             "B54-S1","B55-S1","B56-S1","B57-S1","B58-S1",
                             "B59-S1"
                             };

  string folder[n]={
      "P50-1",      "P50-1",      "P50-1",      "P50-1",      "P50-1",
      "LiLS_batch1_sample1",      "LiLS_batch1_sample2",      "LiLS_batch2_sample1",      "LiLS_batch3_sample1",      "LiLS_batch4_sample1",
      "LiLS_batch5_sample1",      "LiLS_batch6_sample1",      "LiLS_batch7_sample1",      "LiLS_batch8_sample1",      "LiLS_batch9_sample1",
      "LiLS_batch10_sample1",      "LiLS_batch12_sample1",      "LiLS_batch13_sample1",      "LiLS_batch14_sample1",      "LiLS_batch15_sample1",
      "LiLS_batch16_sample1",      "LiLS_batch17_sample1",      "LiLS_batch18_sample1",      "LiLS_batch19_sample1",      "LiLS_batch20_sample1",
      "LiLS_batch21_sample1",      "LiLS_batch22_sample1",      "LiLS_batch23_sample1",      "LiLS_batch24_sample1",      "LiLS_batch25_sample1",
      "LiLS_batch25_sample1_repeat",      "LiLS_batch25_sample2",      "LiLS_drum12_sample1",      "LiLS_batch26_sample1",      "LiLS_batch27_sample1",
      "LiLS_batch28_sample1",      "LiLS_batch29_sample1",      "LiLS_batch30_sample1",      "LiLS_batch31_sample1",      "LiLS_batch32_sample1",
      "LiLS_batch33_sample1",      "LiLS_batch34_sample1",      "LiLS_batch35_sample1",      "LiLS_batch36_sample1",      "LiLS_batch37_sample1",
      "LiLS_batch38_sample1",      "LiLS_batch39_sample1",      "LiLS_batch40_sample1",      "LiLS_batch41_sample1",      "LiLS_batch42_sample1",
      "LiLS_batch43_sample1",      "LiLS_batch44_sample1",      "LiLS_batch45_sample1",      "LiLS_batch47_sample1",      "LiLS_batch48_sample1",
      "LiLS_batch49_sample1",      "LiLS_batch50_sample1",      "LiLS_batch51_sample1",      "LiLS_batch52_sample1",      "LiLS_batch53_sample1",
      "LiLS_batch54_sample1",      "LiLS_batch55_sample1",      "LiLS_batch56_sample1",      "LiLS_batch57_sample1",      "LiLS_batch58_sample1",
      "LiLS_batch59_sample1"
  };
  string AmBeRun[n]={
      "run3568386939",      "run3568565429",      "run3568649727",      "run3569856498",      "run3570538570",
      "run3571322081",      "run3571665388",      "run3571748994",      "run3571835508",      "run3571917776",
      "run3571961506",      "run3572274897",      "run3572299453",      "run3572868578",      "run3572904539",
      "run3573476634",      "run3574807487",      "run3574851408",      "run3574771357",      "run3574683370",
      "run3575292183",      "run3575370827",      "run3575451175",      "run3575893683",      "run3575912813",
      "run3576495873",      "run3577099972",      "run3577110578",      "run3577195912",      "run3577708747",
      "run3577807179",      "run3579524482",      "run3582042442",      "run3577718283",      "run3577727981",
      "run3577739852",      "run3579534033",      "run3579542057",      "run3579552150",      "run3579560793",
      "run3579602567",      "run3579612373",      "run3580125189",      "run3580135031",      "run3580143956",
      "run3581338281",      "run3581942991",      "run3581952181",      "run3581960079",      "run3581978472",
      "run3582020464",      "run3582028068",      "run3582544897",      "run3583756895",      "run3583772667",
      "run3583790193",      "run3583798199",      "run3583836678",      "run3583859630",      "run3584358486",
      "run3584367312",      "run3584375767",      "run3586793005",      "run3586800674",      "run3594657277",
      "run3594666015"
  };


   float FOMA[n], FOMAerr[n];
   float theZ[n], theZerr[n];
   float EperQ[n];
   TH1F* histoZ_P50 = new TH1F("h1000Z_P50","",60,100,130);
   TH1F* histoZ_LiLS = new TH1F("h1000Z_LiLS","",60,100,130);
   TH1F* histoFom_P50 = new TH1F("hFom_P50","",36,1.4-0.025/2,2.3-0.025/2);
   TH1F* histoFom_LiLS = new TH1F("hFom_LiLS","",36,1.4-0.025/2,2.3-0.025/2);
   TH1F* histoEperQ_P50 = new TH1F("hEperQ_P50","",26,22,35);
   TH1F* histoEperQ_LiLS = new TH1F("hEperQ_LiLS","",26,22,35);

   for(int i=0;i<n;i++){
      string filename = folder[i] + "_" + AmBeRun[i] + "_FOM_Z_values.txt";
      fstream fin(filename.c_str(),ios::in);
      while(fin>>total>>fast>>energycal>>fom>>fomerr>>Z>>Zerr){
        if(total==800 && fast==100){
          FOMA[i]=fom; FOMAerr[i]=fomerr;
          theZ[i]=Z; theZerr[i]=Zerr;
          EperQ[i]=energycal;
          break;
        }
      }
      fin.close();
      if(i<=4){
           histoZ_P50->Fill(1000.0*Z);
           histoFom_P50->Fill(fom);
           histoEperQ_P50->Fill(energycal);
      }
      else{
           histoZ_LiLS->Fill(1000.0*Z);
           histoFom_LiLS->Fill(fom);
           histoEperQ_LiLS->Fill(energycal);
      }
   }

/*
   TCanvas *c1 = new TCanvas("cTau1","cTau1",10,10,1640,960);
   c1->SetGrid();
   c1->SetLeftMargin(0.15);
   c1->SetBottomMargin(0.15);

    TH2F* h = new TH2F("h","",n+1,0,n+1,35,22,30);
    TH1F* h0 = new TH1F("h0","h0",n+1,0,n+1);
    h0->SetMarkerStyle(20);
   h->SetStats(0);
   h0->SetStats(0);
   gRandom->SetSeed();
   for (Int_t i=0;i<n;i++) {
      h->Fill(people[i],EperQ[i],EperQ[i]);
      h->SetBinError(i,i,0.000000001);
      h0->SetBinContent(i+1,EperQ[i]);
      h0->SetBinError(i+1,0.000000001);
   }
   h->LabelsDeflate("X");
   //h->LabelsDeflate("Y");
   h->LabelsOption("v");
   h->Draw("same""axis""text");
   h0->Draw("E""same");
   h->SetXTitle("Samples");
   h->GetXaxis()->SetTitleOffset(-0.5);
   h->SetYTitle("Energy calibration constant");
*/
   TCanvas* c2 = new TCanvas();
   histoZ_LiLS->Draw(); histoZ_LiLS->SetLineColor(kRed);
   histoZ_P50->Draw("same"); histoZ_P50->SetLineColor(kBlack);
   TLegend* leg2 = new TLegend(0.6,0.4,0.9,0.6);
   leg2->AddEntry(histoZ_LiLS,"LiLS samples","l");
   leg2->AddEntry(histoZ_P50,"P50 sample 1","l");
   leg2->Draw();

   TCanvas* c3 = new TCanvas();
   histoFom_LiLS->Draw(); histoFom_LiLS->SetLineColor(kRed);
   histoFom_P50->Draw("same"); histoFom_P50->SetLineColor(kBlack);
   histoFom_LiLS->SetXTitle("FOM");
   histoFom_LiLS->SetYTitle("Counts");
   TLegend* leg3 = new TLegend(0.2,0.7,0.5,0.9);
   leg3->AddEntry(histoFom_LiLS,"LiLS samples","l");
   leg3->AddEntry(histoFom_P50,"P50 sample 1","l");
   leg3->Draw();

   TCanvas* c4 = new TCanvas();
   histoEperQ_LiLS->Draw(); histoEperQ_LiLS->SetLineColor(kRed);
   histoEperQ_P50->Draw("same"); histoEperQ_P50->SetLineColor(kBlack);
   histoEperQ_LiLS->SetXTitle("EperQ");
   histoEperQ_LiLS->SetYTitle("Counts");
   TLegend* leg4 = new TLegend(0.2,0.7,0.5,0.9);
   leg4->AddEntry(histoEperQ_LiLS,"LiLS samples","l");
   leg4->AddEntry(histoEperQ_P50,"P50 sample 1","l");
   leg4->Draw();
}
