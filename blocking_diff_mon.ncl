

	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

	begin
	

	inpath      = "/mnt/e/data/"
	
	oupath      = "~/onedrive/EAfront/fig/"

	  latSt     = 40
	  latLt  	= 75
	  lonSt 	= 0
	  lonLt  	= 180
      yrSt  	= 1970
	  yrmd 		= 1992
	  yrLt   	= 2020
      
	  monSt 	= 6
	  month 	= "June"
		  
	  nyr1    	= yrmd - yrSt + 1   ; 1979-1992
	  nyr2    	= yrLt - yrmd		; 1993-2018
    
	  fp       	= addfile(inpath + "monthly_blocking_count_ncep_194801-202006.nc","r")
	 
	  pRaw1 		= fp->blocking_count({yrSt:yrmd},{monSt},:,:)
	  pRaw2 		= fp->blocking_count({yrmd+1:yrLt},{monSt},:,:)
	   	    
	
	  ave1 			:= dim_avg_n_Wrap(pRaw1,0)
	  var1 			:= dim_variance_n_Wrap(pRaw1,0)
	 
	  ave2 			:= dim_avg_n_Wrap(pRaw2,0)
	  var2 			:= dim_variance_n_Wrap(pRaw2,0)
	  
	  prob 		= ttest(ave1,var1,nyr1,ave2,var2,nyr2,False,False)
	  
	  
      dp 		= ave2 - ave1  
	 
		
      copy_VarCoords(pRaw1(0,:,:),dp)
	  copy_VarCoords(pRaw1(0,:,:),prob)
	  
	  
	 printMinMax(dp,0)
	 ;exit
	;==========plot==============
	
	
    wksname 	= oupath+"blocking_count_diff_"+month+"_95"
	
	wks 		= gsn_open_wks("x11",wksname)
	
	gsn_define_colormap(wks,"BlueDarkRed18") 
	plot 		= new(1,graphic)


	res          				= True
	res@gsnDraw 				= False
	res@gsnFrame 				= False
	res@gsnAddCyclic 			= False

	res@cnInfoLabelOn         	= False
	res@cnFillOn              	= True               ; turn on color for contours
	res@cnLinesOn            	= False             ; turn off contour lines
	res@cnLineLabelsOn       	= False              ; turn off contour line labels
    
	; lonSt 		= -100
	; lonLt 		= 180
	res@mpCenterLonF 			= (lonSt+lonLt)/2.
	res@mpMaxLatF 				= latLt
	res@mpMinLatF		 		= latSt
	res@mpMinLonF 				= lonSt
	res@mpMaxLonF 				= lonLt

	res@mpGeophysicalLineColor   	   = "black"
	res@mpGeophysicalLineThicknessF = 1.5

	res@tmXTOn    				= False
	res@tmYROn    				= False
	res@tmXBLabelFontHeightF 	= 0.02
	res@tmYLLabelFontHeightF 	= 0.02
	res@tmXBTickSpacingF  		= 30

	res@gsnSpreadColors 		= True

	res@pmLabelBarOrthogonalPosF  = 0.3
	res@lbLabelFontHeightF       	= 0.015
	res@lbLabelBarOn             	= True
	res@lbBoxLinesOn             	= False
	;res@lbOrientation				= "vertical"

	res@cnLevelSelectionMode 		= "ExplicitLevels"
	res@cnFillColors      			= (/3,4,5,7,9,10,11,12,13,15,17,18/)
	
	res@tiMainFontHeightF 			= 0.02
	res@gsnRightString 			= " "
    
   
 	
	res@cnLevels  			= (/-5,-4,-3,-2,-1,0,1,2,3,4,5/)*0.4

	  
	   res@tiMainString 	= "Diff. Blocking Count in "+month+" (1993-2018 minus 1979-1992)"
	   plot          		= gsn_csm_contour_map_ce(wks,dp,res)
	
	   
	 

	delete(res@mpCenterLonF)
	delete(res@mpMaxLatF)
	delete(res@mpMinLatF)
	delete(res@mpMaxLonF)
	delete(res@mpMinLonF)
	delete(res@mpGeophysicalLineColor)
	delete(res@mpGeophysicalLineThicknessF)
	delete(res@cnFillColors )
	
	
	delete(res@cnLevels)
	res@cnFillOn 						= False
	res@cnLinesOn           			= False              ; turn off contour lines
	res@cnLineLabelsOn       			= False
	res@cnLevels 				       	= (/0.05,0.1/)
	res@cnFillDotSizeF 				= 0.003
	res@cnFillColor 					= "gray40"
	res@gsnLeftString					= " "
	res@tiMainString 					= " "
	shade                  				= True
	shade@gsnShadeFillType 				= "pattern"
	shade@gsnShadeLow 	   				= 17
	;shade@gsnShadeHigh 					= 17
	
	
	
	   plot0         = gsn_csm_contour(wks,prob,res)
	   plot0         = gsn_contour_shade(plot0,0.05,0.1,shade)
	   overlay(plot,plot0)
  
	 
  
	draw(plot)
	frame(wks)
	end


