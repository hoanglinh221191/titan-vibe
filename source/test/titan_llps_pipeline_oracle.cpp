#include "titan_llps.h"
#include "MOL.h"
#include <iostream>
using namespace titan_llps_geometry;
void require(bool x,const char* m){if(!x)throw std::runtime_error(m);}
void close(double a,double b){require(std::isfinite(a)&&std::abs(a-b)<1e-4,"pipeline mismatch");}
int main(){try{
 std::vector<mol_info> metadata;
 for(int i=0;i<4;++i)metadata.emplace_back(i+1,"ALA","BB",std::string(1,char('A'+i)),"C");
 for(const std::string geometry:{"droplet","slab"}) {
  TitanLLPS sync(4,{0,1,2,3},{1,1,1,1},.7,.5,1,false,.5,false,false,0,"scalar","full");
  TitanLLPS async(4,{0,1,2,3},{1,1,1,1},.7,.5,1,false,.5,false,false,0,"auto","full");
  sync.set_sticker_beads(metadata,"BB");async.set_sticker_beads(metadata,"BB");
  sync.configure_phase({geometry,2,1,2,3});async.configure_phase({geometry,2,1,2,3});async.start_async();
  double density_integral_sum=0;
  for(int frame=0;frame<6;++frame) {
    MOL mol;mol.pos.resize(4);mol.Box=box(200.0+frame,200.0+frame,200.0+frame,80.,85.,65.);
    auto cell=Cell::from_box(mol.Box,0,0,0);
    density_integral_sum+=4/cell.volume;
    const Point pos[4]={{0,0,.2},{0,0,-.2},{0,0,2.5},{0,0,7}};
    for(int i=0;i<4;++i){auto d=cell.cartesian({double(frame+i),double(-frame),double(i-frame)});for(int k=0;k<3;++k){auto& v=k==0?mol.pos.x:k==1?mol.pos.y:mol.pos.z;v[i]=10*(pos[i][k]+d[k]);}}
    sync.process_frame(mol,metadata,20,20,20,frame,frame);
    async.process_frame_async(mol,metadata,20,20,20,frame,frame);
  }
  async.flush_async();require(sync.result().frames.size()==6&&async.result().frames.size()==6,"lost frame");
  require(sync.density_frame_count_==6&&async.density_frame_count_==6,"density frame lost");
  double total=0;for(std::size_t i=0;i<sync.vox_pro_.size();++i){total+=sync.vox_pro_[i];close(sync.vox_pro_[i],async.vox_pro_[i]);}
  close(total/sync.vox_pro_.size(),density_integral_sum);
  for(int i=0;i<6;++i){auto a=sync.result().frames[i],b=async.result().frames[i];require(a.gcc_size==2&&a.phase.status=="ok"&&a.phase.dense_count==2&&a.phase.dilute_count==1,"wrong phase pipeline counts");close(a.gcc_rg_nm,.2);close(a.density_dense,b.density_dense);close(a.density_dilute,b.density_dilute);close(a.gcc_rg_nm,b.gcc_rg_nm);}
 }
 TitanLLPS winding(4,{0,1,2,3},{1,1,1,1},.6,.5,1,false,.5,false,false,0,"auto","off");
 winding.set_sticker_beads(metadata,"BB");winding.configure_phase({"droplet",2,.1,.2,.3});
 MOL ring;ring.pos.resize(4);ring.Box=box(20.,20.,20.,90.,90.,90.);for(int i=0;i<4;++i)ring.pos.x[i]=5*i;
 winding.process_frame(ring,metadata,2,2,2,0,0);
 require(std::isnan(winding.result().frames[0].gcc_rg_nm)&&winding.result().frames[0].phase.status=="gcc_embedding_ambiguous","winding graph must abstain");
 TitanLLPS prna(4,{0,1,2,3},{1,1,1,1},.7,.5,1,true,.5,true,false,0,"auto","off");
 prna.set_sticker_beads(metadata,"BB");std::vector<mol_info> rna_info{mol_info(1,"P","P","R","RNA")};prna.set_rna_metadata(rna_info);
 MOL pro,rna;pro.pos.resize(4);rna.pos.resize(1);pro.Box=box(200.,200.,200.,80.,85.,65.);rna.Box=pro.Box;
 pro.pos.z[0]=-2;pro.pos.z[1]=2;pro.pos.z[2]=25;pro.pos.z[3]=70;
 auto cell=Cell::from_box(pro.Box,0,0,0);auto shift=cell.cartesian({2,-3,4});rna.pos.x[0]=shift[0]*10;rna.pos.y[0]=shift[1]*10;rna.pos.z[0]=shift[2]*10+3;
 prna.process_frame(pro,metadata,20,20,20,0,0,&rna,&rna_info);
 require(prna.protein_rna_pair_keys_buf_.size()==2,"triclinic protein/RNA contacts");
 TitanLLPS bad(4,{0,1,2,3},{1,1,1,1});bad.set_sticker_beads(metadata,"BB");
 bad.start_async();ring.pos.x[0]=std::numeric_limits<double>::quiet_NaN();
 bool caught=false;try{bad.process_frame_async(ring,metadata,2,2,2,0,0);bad.flush_async();}catch(const std::invalid_argument&){caught=true;}
 require(caught,"async worker exception lost");
 std::cout<<"PASS LLPS pipeline: droplet/slab synchronous/async, changing triclinic cells, image shifts, winding rejection, worker exception propagation\n";
 return 0;
}catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
