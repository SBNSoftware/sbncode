namespace single_photon
{

    //ask YJ
    int SinglePhoton::setTPCGeom(){

        m_tpc_active_x_low  =  0.0;
        m_tpc_active_x_high =  256.35;
        m_tpc_active_y_low  =  -116.5;
        m_tpc_active_y_high =  116.5;
        m_tpc_active_z_low  =  0.0;
        m_tpc_active_z_high =  1036.8;

        m_SCB_YX_TOP_y1_array =  116.;
        m_SCB_YX_TOP_x1_array =  {0., 150.00, 132.56, 122.86, 119.46, 114.22, 110.90, 115.85, 113.48, 126.36, 144.21};
        m_SCB_YX_TOP_y2_array =  {0., 110.00, 108.14, 106.77, 105.30, 103.40, 102.18, 101.76, 102.27, 102.75, 105.10};
        m_SCB_YX_TOP_x2_array =  256.;

        m_SCB_YX_BOT_y1_array = -115.;
        m_SCB_YX_BOT_x1_array = {0., 115.71, 98.05, 92.42, 91.14, 92.25, 85.38, 78.19, 74.46, 78.86, 108.90};
        m_SCB_YX_BOT_y2_array = {0., -101.72, -99.46, -99.51, -100.43, -99.55, -98.56, -98.00, -98.30, -99.32, -104.20};
        m_SCB_YX_BOT_x2_array = 256.;

        m_SCB_ZX_Up_z1_array = 0.;
        m_SCB_ZX_Up_x1_array = 120.;
        m_SCB_ZX_Up_z2_array = 11.;
        m_SCB_ZX_Up_x2_array = 256.;//upstream

        m_SCB_ZX_Dw_z1_array     = 1037.;
        m_SCB_ZX_Dw_x1_array = {0., 120.00, 115.24, 108.50, 110.67, 120.90, 126.43, 140.51, 157.15, 120.00, 120.00};
        m_SCB_ZX_Dw_z2_array = {0., 1029.00, 1029.12, 1027.21, 1026.01, 1024.91, 1025.27, 1025.32, 1027.61, 1026.00, 1026.00};
        m_SCB_ZX_Dw_x2_array     = 256.;//downstream

        return 0;
    }


    int SinglePhoton::isInTPCActive(std::vector<double> & vec){
        return isInTPCActive(0.0,vec);
    }

    /* determine if point vec is inside TPC active volume, returns 1 - in TPC active, 0 - out of TPC active */
    int SinglePhoton::isInTPCActive(double cut,std::vector<double> & vec){
        bool is_x = (vec[0] > m_tpc_active_x_low+cut && vec[0]< m_tpc_active_x_high-cut );
        bool is_y = (vec[1] > m_tpc_active_y_low+cut && vec[1]< m_tpc_active_y_high-cut);
        bool is_z = (vec[2] > m_tpc_active_z_low+cut && vec[2]<m_tpc_active_z_high-cut );
        return ( (is_x && is_y && is_z) ? 1 : 0  );

    }

    /* returns minimum distance to the TPC active boundary; returns -999 if the point is not in TPC active volume */
    double SinglePhoton::distToTPCActive(std::vector<double>&vec){
        if(isInTPCActive(vec)==0) return -999;
        double min_x = std::min( fabs(vec[0] - m_tpc_active_x_low) ,  fabs(vec[0] - m_tpc_active_x_high));
        double min_y = std::min( fabs(vec[1] - m_tpc_active_y_low) ,  fabs(vec[1] - m_tpc_active_y_high));
        double min_z = std::min( fabs(vec[2] - m_tpc_active_z_low) ,  fabs(vec[2] - m_tpc_active_z_high));

        return ( (min_x<min_y) ? std::min(min_x,min_z)  : std::min(min_y,min_z)    );
    }


    /* Guanqun: load the SCB polygon in XY plane at various Z position ? */
    bool SinglePhoton::loadSCB_YX(std::vector<TGeoPolygon*>& zpolygons){
        //TGeoManager *geom = new TGeoManager("save scb", "save scb");
        //cout << "size of " << zpolygons.size() << endl;

        //DEPRECIATED SHOULD NOT BE USED!!!
        std::cout<<"ERROR ERROR ERROR DEPRECIATED DEPRECIATED. "<<std::endl;

        double tbi = -10000.; //means "to be initialized"                                                                                                                                
        //torso
        double ptX[6] = {0., tbi, m_SCB_YX_TOP_x2_array, m_SCB_YX_BOT_x2_array, tbi, 0.};
        double ptY[6] = {m_SCB_YX_TOP_y1_array, m_SCB_YX_TOP_y1_array, tbi, tbi, m_SCB_YX_BOT_y1_array, m_SCB_YX_BOT_y1_array};

        TGeoPolygon *polyXY = new TGeoPolygon(6);

        for (Int_t z_idx_YX=0; z_idx_YX<10; z_idx_YX++) {
            ptX[1] = m_SCB_YX_TOP_x1_array[z_idx_YX+1];
            ptX[4] = m_SCB_YX_BOT_x1_array[z_idx_YX+1];
            ptY[2] = m_SCB_YX_TOP_y2_array[z_idx_YX+1];
            ptY[3] = m_SCB_YX_BOT_y2_array[z_idx_YX+1];
            polyXY->SetXY(ptX,ptY);
            polyXY->FinishPolygon();
            zpolygons.push_back(polyXY); // Guanqun: not sure pushing pack pointers works as expect?? 
        }

        //cout << "size of " <<zpolygons.size() << endl;
        return 0;
    }


    int SinglePhoton::isInSCB(std::vector<double> & vec){
        return isInSCB(0.0,vec);
    }

    int SinglePhoton::distToSCB(double & dist, std::vector<double> &vec){
        //this one returns the distance to the boundary
        bool ans = false;
        double dist_yx = 999999;

        int iseg = 0;
        if (!isInTPCActive(vec)){
            dist=-1;
            return 0; // is it in active volume?
        }
        double cut = 0.0;

        TVector3 pt(&vec[0]);

        Int_t z_idx = (Int_t)pt.Z()/100;  // Int_t: signed integer 4 bytes

        Int_t z_idx_YX = z_idx;//YX-view effective z index, it is only different for z > 10m area, where we want to appliy 9m<z<10m YX boundary, still need to keep the original z_idx bc it's needed in ZX-view
        if (z_idx_YX==10) z_idx_YX-=1;
        double tbi = -10000.; //means "to be initialized"


        double ptX[6] = {0.+cut, tbi, m_SCB_YX_TOP_x2_array-cut, m_SCB_YX_BOT_x2_array-cut, tbi, 0.+cut};
        double ptY[6] = {m_SCB_YX_TOP_y1_array-cut, m_SCB_YX_TOP_y1_array-cut, tbi, tbi, m_SCB_YX_BOT_y1_array+cut, m_SCB_YX_BOT_y1_array+cut};

        TGeoPolygon *polyXY = new TGeoPolygon(6);

        ptX[1] = m_SCB_YX_TOP_x1_array[z_idx_YX+1];
        ptX[4] = m_SCB_YX_BOT_x1_array[z_idx_YX+1];
        ptY[2] = m_SCB_YX_TOP_y2_array[z_idx_YX+1];
        ptY[3] = m_SCB_YX_BOT_y2_array[z_idx_YX+1];

        polyXY->SetXY(ptX,ptY);
        polyXY->FinishPolygon();
        double testpt[2] = {pt.X(), pt.Y()};

        //cout << "is testpt ("<< pt.X()<<", "<<pt.Y()<<") contrained? "<<  polyXY->Contains(testpt)<<endl;
        //cout << "area ? " << polyXY->Area()<< endl;

        //polyXY->Draw();    

        Bool_t XY_contain = polyXY->Contains(testpt);
        dist_yx = polyXY->Safety(testpt, iseg); // Compute minimum distance from testpt to any segment.

        if(0<z_idx && z_idx<10){
            double up_z = pt.Z()-m_tpc_active_z_low; // gonna bet, if it's middle enough to be up_z or down_z is smaller than the safefy in this z regime (1m,10m), it is safe to set up_z = z-0, down_z=1036.8-z 
            double down_z = m_tpc_active_z_high-pt.Z();
            double min_z =  std::min(up_z,down_z);

            XY_contain ?  dist = std::min(dist_yx,min_z) : dist = -1  ;

            delete polyXY;
            return  (XY_contain ? 1 : 0 );
        }

        //up or down
        double top_y = m_tpc_active_y_high-pt.Y();
        double bottom_y = pt.Y()+m_tpc_active_y_high;
        double min_y = std::min(top_y, bottom_y);


        // if z_idx==0 or z_idx==10, they need xz view,  
        /// ZX view has Y dependence: Y sub-range from -116 to 116cm per 24cm
        Int_t y_idx = (pt.Y()+116.)/24;
        if (pt.Y()<-116. && pt.Y()>-116.5) y_idx = 0; //just the 0.5cm 
        if(y_idx<0 || y_idx>9) {
            dist = -1; 
            delete polyXY;  
            return 0;
        }

        Bool_t ZX_contain = false;
        double arbit_out = 55555;

        double d_dist = -1;

        //upstream
        if(z_idx==0){
            double ptX_Up[5] = {0.+cut, m_SCB_ZX_Up_x1_array, m_SCB_ZX_Up_x2_array-cut, m_SCB_ZX_Up_x2_array-cut, 0+cut};
            double ptZ_Up[5] = {0.+cut,0.+cut,m_SCB_ZX_Up_z2_array, arbit_out, arbit_out};

            TGeoPolygon *polyXZ_Up = new TGeoPolygon(5);
            polyXZ_Up->SetXY(ptX_Up,ptZ_Up);
            polyXZ_Up->FinishPolygon();

            double testpt_Up[2] = {pt.X(), pt.Z()};
            ZX_contain = polyXZ_Up->Contains(testpt_Up);
            d_dist = polyXZ_Up->Safety(testpt_Up,iseg);
            delete polyXZ_Up;
   
        }
        if (z_idx==10){
            double ptX_Dw[5] = {0.+cut,m_SCB_ZX_Dw_x2_array-cut, m_SCB_ZX_Dw_x2_array-cut, tbi, 0.+cut};
            double ptZ_Dw[5] = {-arbit_out,-arbit_out,tbi,m_SCB_ZX_Dw_z1_array-cut, m_SCB_ZX_Dw_z1_array-cut};

            ptX_Dw[3] = m_SCB_ZX_Dw_x1_array[y_idx+1];
            ptZ_Dw[2] = m_SCB_ZX_Dw_z2_array[y_idx+1];

            TGeoPolygon *polyXZ_Dw = new TGeoPolygon(5);
            polyXZ_Dw->SetXY(ptX_Dw,ptZ_Dw);
            polyXZ_Dw->FinishPolygon();

            double testpt_Dw[2] = {pt.X(), pt.Z()};
            ZX_contain = polyXZ_Dw->Contains(testpt_Dw);
            d_dist = polyXZ_Dw->Safety(testpt_Dw,iseg);
            delete polyXZ_Dw;
        }

        delete polyXY;

        ans = XY_contain && ZX_contain;
        ans ? dist = std::min(d_dist,min_y) :  dist=-1;

        return (ans ? 1 : 0);
    }


    int SinglePhoton::isInSCB(double cut, std::vector<double> &vec){

        if (!isInTPCActive(cut,vec)) return 0; // is it in active volume?
        TVector3 pt(&vec[0]);

        Int_t z_idx = (Int_t)pt.Z()/100;

        Int_t z_idx_YX = z_idx;//YX-view effective z index, it is only different for z > 10m area, where we want to appliy 9m<z<10m YX boundary, still need to keep the original z_idx bc it's needed in ZX-view
        if (z_idx_YX==10) z_idx_YX-=1;
        double tbi = -10000.; //means "to be initialized"

        double ptX[6] = {0.+cut, tbi, m_SCB_YX_TOP_x2_array-cut, m_SCB_YX_BOT_x2_array-cut, tbi, 0.+cut};
        double ptY[6] = {m_SCB_YX_TOP_y1_array-cut, m_SCB_YX_TOP_y1_array-cut, tbi, tbi, m_SCB_YX_BOT_y1_array+cut, m_SCB_YX_BOT_y1_array+cut};

        TGeoPolygon *polyXY = new TGeoPolygon(6);

        ptX[1] = m_SCB_YX_TOP_x1_array[z_idx_YX+1];
        ptX[4] = m_SCB_YX_BOT_x1_array[z_idx_YX+1];
        ptY[2] = m_SCB_YX_TOP_y2_array[z_idx_YX+1];
        ptY[3] = m_SCB_YX_BOT_y2_array[z_idx_YX+1];

        polyXY->SetXY(ptX,ptY);
        polyXY->FinishPolygon();
        double testpt[2] = {pt.X(), pt.Y()};

        //cout << "is testpt ("<< pt.X()<<", "<<pt.Y()<<") contrained? "<<  polyXY->Contains(testpt)<<endl;
        //cout << "area ? " << polyXY->Area()<< endl;

        //polyXY->Draw();    

        Bool_t XY_contain = polyXY->Contains(testpt);

        if(0<z_idx && z_idx<10){
            delete polyXY;
            return XY_contain;
        }

        // if z_idx==0 or z_idx==10, they need xz view,  

        /// ZX view has Y dependence: Y sub-range from -116 to 116cm per 24cm
        Int_t y_idx = (pt.Y()+116.)/24;

        if (pt.Y()<-116. && pt.Y()>-116.5) y_idx = 0; //just the 0.5cm 

        if(y_idx<0 || y_idx>9) {
            delete polyXY;  
            return 0;
        }

        Bool_t ZX_contain = false;

        if(z_idx==0){
            double ptX_Up[5] = {0.+cut, m_SCB_ZX_Up_x1_array, m_SCB_ZX_Up_x2_array-cut, m_SCB_ZX_Up_x2_array-cut, 0+cut};
            double ptZ_Up[5] = {0.+cut,0.+cut,m_SCB_ZX_Up_z2_array, 100., 100.};

            TGeoPolygon *polyXZ_Up = new TGeoPolygon(5);
            polyXZ_Up->SetXY(ptX_Up,ptZ_Up);
            polyXZ_Up->FinishPolygon();

            double testpt_Up[2] = {pt.X(), pt.Z()};
            ZX_contain = polyXZ_Up->Contains(testpt_Up);
            delete polyXZ_Up;
        }

        else if (z_idx==10){
            double ptX_Dw[5] = {0.+cut,m_SCB_ZX_Dw_x2_array-cut, m_SCB_ZX_Dw_x2_array-cut, tbi, 0.+cut};
            double ptZ_Dw[5] = {1000.,1000.,tbi,m_SCB_ZX_Dw_z1_array-cut, m_SCB_ZX_Dw_z1_array-cut};

            ptX_Dw[3] = m_SCB_ZX_Dw_x1_array[y_idx+1];
            ptZ_Dw[2] = m_SCB_ZX_Dw_z2_array[y_idx+1];

            TGeoPolygon *polyXZ_Dw = new TGeoPolygon(5);
            polyXZ_Dw->SetXY(ptX_Dw,ptZ_Dw);
            polyXZ_Dw->FinishPolygon();

            double testpt_Dw[2] = {pt.X(), pt.Z()};
            ZX_contain = polyXZ_Dw->Contains(testpt_Dw);
            delete polyXZ_Dw;
        }

        delete polyXY;
        return (XY_contain && ZX_contain);

    }

}
