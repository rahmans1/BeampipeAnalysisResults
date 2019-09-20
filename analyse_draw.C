#include "remolltypes.hh"
using namespace ROOT;

int analyse(TString conf, Int_t nFiles, TString particle, TString energy){

TChain T("T");

T.Add(Form("/scratch/j/jmammei/rahmans/scratch/targetStudy/shortened_newcol12_uv_JLAB_75/beam_20mil_%s/beam_%d.root",conf.Data(), nFiles));

Int_t nEvents=10*T.GetEntries();
//Int_t nEventsPerCore=(int) nEvents/nCore;
//Int_t start= (int) (nJob-1)*nEventsPerCore;
//Int_t end= (int) nJob*nEventsPerCore;
//Double_t weight= (1.0/nEventsPerCore)*(1.0/(1.602*1e-13));
std::cout<< "Analyzing  "<<nEvents<<" from " << conf.Data()<<"_"<<nFiles << std::endl;

TFile f(Form("/scratch/j/jmammei/rahmans/root/beampipeStudy/new/%s_%d_%s_%s_%i.root",conf.Data(),nEvents,particle.Data(),energy.Data(),nFiles), "RECREATE");


/*
TH2D rz("rz", "hit.r vs hit.z, bpipe [Hz/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
TH2D rz_pdep("rz_pdep", "(hit.r vs hit.z)*hit.edep, bpipe [W/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
TH2D rz_pin("rz_pin", "(hit.r vs hit.z)*hit.p, bpipe [W/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
TH2D vrvz("vrvz", "hit.vr vs hit.vz, bpipe [Hz/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
TH2D vrvz_pdep("vrvz_pdep", "(hit.vr vs hit.vz)*hit.edep, bpipe [W/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
TH2D vrvz_pin("vrvz_pin", "(hit.vr vs hit.vz)*hit.p, bpipe [W/uA/mm^2]", 29500, 0, 29500, 250, 0, 250);
*/

TH2D rvz("rvz","(hit.r vs hit.vz), det [Hz/uA/mm^2]", 29500, 0, 29500, 1900, 0, 1900);
TH2D rvz_p("rvz_p","(hit.r vs hit.vz)*hit.p, det [W/uA/mm^2]", 29500, 0, 29500, 1900, 0, 1900);

TH2D rvz_rcut("rvz_rcut", "(hit.r vs hit.vz), det and 50cm<=hit.r<=190cm [Hz/uA/mm^2]", 29500, 0, 29500, 1400, 500, 1900);
TH2D rvz_rcut_p("rvz_rcut_p", "(hit.r vs hit.vz)*hit.p, det and 50cm<=hit.r<=190cm [W/uA/mm^2]", 29500, 0, 29500, 1400, 500, 1900);

TH2D vrvz_rcut("vrvz_rcut", "(hit.vr vs hit.vz), det and 50cm<=hit.r<=190cm  [Hz/uA/mm^2]", 29500, 0, 29500, 800, 0, 800);
TH2D vrvz_rcut_p("vrvz_rcut_p", "(hit.vr vs hit.vz)*hit.p, det and 50cm<=hit.r<=190cm [W/uA/mm^2]", 29500, 0, 29500, 800, 0, 800);

TH2D xy_rcut("xy_rcut", "(hit.x vs hit.y), det and 50cm<=hit.r<=190cm [Hz/uA/mm^2]", 3800, -1900,1900, 3800,-1900,1900);
TH2D xy_rcut_p("xy_rcut_p", "(hit.x vs hit.y)*hit.p, det and 50cm<=hit.r<=190cm [W/uA/mm^2]", 3800, -1900,1900, 3800,-1900,1900);


std::map<TString, Float_t> hit_e_max={{"all", 11000}, {"lowene",10}, {"midene",100}, {"highene", 11000}};
std::map<TString, Float_t> hit_e_min={{"all", 0}, {"lowene",0}, {"midene",10}, {"highene", 100}};
std::map<TString, Float_t> hit_e_bin={{"all", 10}, {"lowene",0.1}, {"midene",1}, {"highene", 10}};
TH1D ke_rcut("ke_rcut", Form("hit.p, det and 50cm<=hit.r<=190cm [Hz/uA/%3.3fMeV]",hit_e_bin[energy]), (Int_t) (hit_e_max[energy]-hit_e_min[energy])/hit_e_bin[energy], hit_e_min[energy], hit_e_max[energy]);
TH1D r_rcut("r_rcut", "hit.r, det and 50cm<=hit.r<=190cm [Hz/uA/mm]", 1400,500,1900); 
TH1D r("r","hit.r, det [Hz/uA/mm]",1900,0,1900);


TCut cutoff="hit.p>=1";
std::map<TString, TCut> hit_energy={{"all", "hit.p>=1"}, {"lowene", "hit.p>=1 && hit.p<10"}, {"midene", "hit.p>=10 && hit.p<100"},{"highene", "hit.p>=100"}};
std::map<TString, TCut> hit_pid={{"all", "1==1"}, {"photon", "hit.pid==22"}, {"electron","hit.pid==11 || hit.pid==-11"}, {"neutron", "hit.pid==2112"}};
TCut plaindet= "hit.det==32";
TCut bpipe= "hit.det>=50 && hit.det<=54";
TCut radial="hit.r>=500 && hit.r<=1900";
TCut pe= hit_energy[energy]&& hit_pid[particle];
TCut weight= Form("1/(1.602*1e-13)*1/%i", nEvents);
TCut einw= "hit.p*1.602*1e-13";
TCut edepw="hit.edep*1.602*1e-13";

T.Draw("hit.r>> r", (cutoff&& plaindet && pe)*weight, "goff");
T.Draw("hit.r>> r_rcut", (cutoff && plaindet && pe && radial)*weight, "goff");
T.Draw("hit.p>> ke_rcut", (cutoff&& plaindet && pe && radial)*weight, "goff");

T.Draw("hit.x:hit.y>> xy_rcut", (cutoff && plaindet && pe && radial)*weight, "goff");
T.Draw("hit.x:hit.y>> xy_rcut_p", (cutoff && plaindet && pe && radial)*weight*einw, "goff");



T.Draw("sqrt(hit.vx*hit.vx+hit.vy*hit.vy):hit.vz>> vrvz_rcut", (cutoff && plaindet && pe&&  radial)*weight, "goff");
T.Draw("sqrt(hit.vx*hit.vx+hit.vy*hit.vy):hit.vz>> vrvz_rcut_p", (cutoff && plaindet && pe&&  radial)*weight*einw, "goff");

T.Draw("hit.r:hit.vz>> rvz", (cutoff && plaindet && pe)*weight, "goff");
T.Draw("hit.r:hit.vz>> rvz_p", (cutoff && plaindet && pe)*weight*einw, "goff");

T.Draw("hit.r:hit.vz>> rvz_rcut", (cutoff && plaindet && pe && radial)*weight, "goff");
T.Draw("hit.r:hit.vz>> rvz_rcut_p",(cutoff && plaindet && pe && radial)*weight*einw, "goff");


/*
T.Draw("sqrt(hit.vx*hit.vx+hit.vy*hit.vy):hit.vz>> vrvz", (cutoff && bpipe && pe)*weight, "goff");
T.Draw("sqrt(hit.vx*hit.vx+hit.vy*hit.vy):hit.vz>> vrvz_pin", (cutoff && bpipe && pe)*weight*einw, "goff");
T.Draw("sqrt(hit.vx*hit.vx+hit.vy*hit.vy):hit.vz>> vrvz_pdep", (cutoff && bpipe && pe)*weight*edepw, "goff");

T.Draw("hit.r:hit.z>> rz", (cutoff && bpipe && pe)*weight, "goff");
T.Draw("hit.r:hit.z>> rz_pin", (cutoff && bpipe && pe)*weight*einw, "goff");
T.Draw("hit.r:hit.z>> rz_pdep", (cutoff && bpipe && pe)*weight*edepw, "goff");
*/


/*
std::vector<remollGenericDetectorHit_t>  *fHit=0;
T.SetBranchAddress("hit", &fHit);

for (size_t j=start;j< end;j++){
	if (j%10000==0) {std::cout<<j<< std::endl;}
	T.GetEntry(j);
        for (size_t i=0;i<fHit->size();i++){
                remollGenericDetectorHit_t hit=fHit->at(i);
                Bool_t hit_beampipe =  hit.det>=50 && hit.det<=54;// && hit.det<=54;
                Bool_t hit_planedet = hit.det==32 ;
                Bool_t hit_radial = hit.r >=500 && hit.r<=1900 ;
                Bool_t hit_cutoff = hit.p<1;    // Cut off all particles with energy less than 1 MeV
                
		std::map<TString, Bool_t> hit_pid={{"all",1},{"electron",hit.pid==11||hit.pid==-11},{"photon", hit.pid==22}, {"neutron",hit.pid==2112}};
                std::map<TString, Bool_t> hit_ene={{"all",1},{"lowene", hit.p>=1 && hit.p<10}, {"midene",  hit.p>=10 && hit.p<100}, {"highene", hit.p>=100}};
                if (hit_cutoff || !hit_pid[particle] || !hit_ene[energy]){
                continue;               
		}

                if (hit_planedet){
                         rvz.Fill(hit.vz, hit.r, weight);
			 rvz_p.Fill(hit.vz, hit.r, 1.602*1e-13*weight*hit.p);
                         r.Fill(hit.r,weight);
 	               if(hit_radial){
                         rvz_rcut.Fill(hit.vz, hit.r, weight);
                         rvz_rcut_p.Fill(hit.vz, hit.r, 1.602*1e-13*weight*hit.p);
                         vrvz_rcut.Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),weight); 
                         vrvz_rcut_p.Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),1.602*1e-13*weight*hit.p); 
                         ke_rcut.Fill(hit.p, weight);
                         r_rcut.Fill(hit.r, weight);
                         xy_rcut.Fill(hit.x,hit.y, weight);
                         xy_rcut_p.Fill(hit.x,hit.y, 1.602*1e-13*weight*hit.p);
        	       }
		} 
                if (hit_beampipe){
                         rz.Fill(hit.z,hit.r, weight);
                         rz_pdep.Fill(hit.z, hit.r, 1.602*1e-13*weight*hit.edep);
                         rz_pin.Fill(hit.z,hit.r, 1.602*1e-13*weight*hit.p);
                         vrvz.Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy), weight);
                         vrvz_pdep.Fill(hit.vz, sqrt(hit.vx*hit.vx+hit.vy*hit.vy), 1.602*1e-13*weight*hit.edep);
                         vrvz_pin.Fill(hit.vz, sqrt(hit.vx*hit.vx+hit.vy*hit.vy), 1.602*1e-13*weight*hit.p);
		}

       }
}
*/

rvz.Write("",TObject::kOverwrite);
rvz_p.Write("", TObject::kOverwrite);
rvz_rcut.Write("", TObject::kOverwrite);
rvz_rcut_p.Write("", TObject::kOverwrite);
r.Write("",TObject::kOverwrite);
/*
vrvz.Write("", TObject::kOverwrite);
vrvz_pdep.Write("", TObject::kOverwrite);
vrvz_pin.Write("", TObject::kOverwrite);

rz.Write("", TObject::kOverwrite);
rz_pdep.Write("", TObject::kOverwrite);
rz_pin.Write("", TObject::kOverwrite);
*/

vrvz_rcut.Write("", TObject::kOverwrite);
vrvz_rcut_p.Write("", TObject::kOverwrite);
ke_rcut.Write("", TObject::kOverwrite);
r_rcut.Write("",TObject::kOverwrite);
xy_rcut.Write("",TObject::kOverwrite);
xy_rcut_p.Write("",TObject::kOverwrite);
return 0;
}
