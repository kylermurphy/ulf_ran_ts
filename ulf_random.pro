function neg_ticks, axis, index, value
  return, string(-1*value,format='(F4.1)')
end


pro ulf_random, sdate, edate, $
  l_min=l_min, $
  l_max=l_max, $
  l_bin=l_bin, $
  mlt_min=mlt_min, $
  mlt_max=mlt_max, $
  mlt_bin=mlt_bin, $
  r_seed = r_seed, $ ; random seed so results are reproducible
  d_dir = d_dir, $
  o_dir = o_dir, $
  w_dat = w_dat, $
  verbose=verbose, $
  movie=movie
  
  
  if keyword_set(l_min) then l_min=l_min else l_min=3.
  if keyword_set(l_max) then l_max=l_max else l_max=7.
  if keyword_set(l_bin) then l_bin=l_bin else l_bin=0.5
  if keyword_set(mlt_min) then mlt_min=mlt_min else mlt_min=0
  if keyword_set(mlt_max) then mlt_max=mlt_max else mlt_max=23
  if keyword_set(mlt_bin) then mlt_bin=mlt_bin else mlt_bin=1
  if keyword_set(r_seed) then r_seed=r_seed else r_seed=17
  if keyword_set(d_dir) then d_dir=d_dir else d_dir='D:\data\magnetometer\psd\'
  if keyword_set(o_dir) then o_dir=o_dir else o_dir='D:\out_put\ulf_ran_ts\'
  if keyword_set(w_dat) then w_dat=1 else w_dat=0
  if keyword_set(verbose) then verbose=verbose else verbose=0
  if keyword_set(movie) then movie=1 else movie=0
  
  ; create l grid
  l_grid = findgen((l_max-l_min)/l_bin + 1)*l_bin+l_min
  ; create centered mlt grid
  m_grid = findgen((mlt_max-mlt_min)/mlt_bin + 1)*mlt_bin+mlt_min+mlt_bin/2.
  
  ; convert l/mlt grid to x/y
  theta = m_grid*(360/24.)*!pi/180.
  x_grid = l_grid # cos(theta)
  y_grid = l_grid # sin(theta)
  
  ; transform to 1-D array
  x_grid = x_grid[*]
  y_grid = y_grid[*]
  
  ; create an hourly time series 
  ; for readin in the PSD files
  sd = time_double(sdate)
  ed = time_double(edate)
  t_arr = sd+dindgen((ed-sd)/3600.)*3600.+1
  t_arr = time_string(t_arr, tformat='YYYY-MM-DD/hh')
  t_arr = time_double(t_arr)
  
  y_i = median(long(time_string(t_arr, tformat='YYYY')))
  ;load station values
  stn_vals = stn_dat(y_i)
  ;load electric field mapping values
  e_map = gb2ee( )
  npts = n_elements(e_map.freq)+1
  
  ;list to hold the final timeseries
  tf = list()
  
  ;get omni data 
  o1h = get_omni_1hr(date_min=sd,date_max=ed)
  x_tk = time_ticks([sd,ed], offset)
  x_tk.xtickv=x_tk.xtickv+offset
  x_tk.xrange=x_tk.xrange+offset
  
  ;movie output
  if keyword_set(movie) then begin
    video_file = time_string(sd,tformat='YYYYMMDD')+'_'+time_string(ed,tformat='YYYYMMDD')
    video_file = 'ULF_ran_ts_'+video_file
    video_file = video_file+'.mp4'
    framerate=5
    ; Open the video recorder.
    video = IDLffVideoWrite(o_dir+video_file, Format='mp4')
    ; Configure the video output for the PNG files you
    ; plan to add.
    stream = video.AddVideoStream(500, 1000, framerate)
  endif
  
  for i=0L, t_arr.length-1 do begin
    fs = d_dir+time_string(t_arr[i],tformat='YYYY/MM/DD/')
    fs = fs+time_string(t_arr[i],tformat='YYYYMMDD')+'*_psd.txt.gz'
    fs = file_search(fs, count=fc)
    
    stn = strsplit(fs,time_string(t_arr[i],tformat='YYYYMMDD'),/regex,/extract)
    stn = stn.ToArray()
    stn = reform(stn[*,-1])
    stn = strsplit(stn,'_',/extract)
    stn = stn.ToArray()
    stn = reform(stn[*,0])

    lsh  = fltarr(stn.length)
    mlt  = fltarr(stn.length)
    b_spec = fltarr(stn.length, npts)
    e_spec = b_spec
    
    h_i = long(time_string(t_arr[i], tformat='hh'))
    r_c = 0L
    for j=0L, stn.length-1 do begin 
      ; get station position
      s_p = where(stn_vals.stn eq stn[j], s_c)
      if s_c ne 1 then continue
      m_val = h_i-stn_vals.mlt_midnight[s_p]
      if m_val lt 0 then m_val = m_val+24
      if m_val gt 24 then m_val = m_val-24
      
      r_s = read_ulf_psd_spec(stn[j], time_string(t_arr[i],tformat='YYYY-MM-DD/hh:00:00'), fmin=-1, fmax=6*1000.)
      
      b_d = where(finite(r_s.psd[h_i,0:npts-1]) ne 1, b_c)
      if b_c gt 0 then continue
      b_t = total(r_s.psd[h_i,0:npts-1])
      if b_t eq 0 then continue
      
      ; convert b psd to nT^2/mHz
      b_spec[r_c,*] = r_s.psd[h_i,0:npts-1]/1000.
      lsh[r_c] = stn_vals.lshell[s_p]
      mlt[r_c] = m_val
      
      ;mapping position
      e_m = min(abs(e_map.l_v-lsh[r_c]))
      e_p = !C
      ;mapping e field
      e1 = reform(b_spec[r_c,*]) ;mapping across all frequencies, backout coeffecient
      e1[0] = 0
      e2 = e1[0:n_elements(e_map.freq)] ;mapping across only frequencies in the file
      
      e1[1:-1] = e1[1:-1]*e_map.rgf2_c[e_p]*(r_s.f[1:-1]^2.)
      e2[1:n_elements(e_map.freq)] = reform(reform(e2[1:n_elements(e_map.freq)]) * reform(e_map.g2e[e_p,*]))
      e2 = reform(e2)
      
      ;check difference between the two types
      e_diff = total(e1[0:e2.length-1]-e2[0:-1])

      e_spec[r_c,*] = e1
      
      
      if verbose gt 1 then print, strtrim(h_i,2)+' '+stn[j]+' '+strtrim(m_val,2)+' '+strtrim(lsh[r_c],2)+' '+strtrim(e_diff,2)
      r_c++
    endfor

    ; only use frequencies upto 10 mHz (the max from the mapping)
    ; this is what Louis recommends
    ; the max frequency from the e_map is 10 mHz
    ; convert the spectra back to [unit]/Hz
    b_spec = b_spec[0:r_c-1,0:n_elements(e_map.freq)]*1000.
    e_spec = e_spec[0:r_c-1,0:n_elements(e_map.freq)]*1000.
    lsh = lsh[0:r_c-1]
    mlt = mlt[0:r_c-1]
    
    ; get nyquist frequency and time step
    nyq = max(e_map.freq)/1000.
    d_t = 1/(2.*nyq)
    
    t_d = mlt*(360/24.)*!pi/180.
    
    x_d = lsh*cos(t_d)
    y_d = lsh*sin(t_d)
    
    z_d = alog10(e_spec[*,1])
    
    ;for loop here for interpolating
    z_i = griddata(x_d,y_d,z_d, method='Kriging', xout=x_grid,yout=y_grid) 
    
    z_spec = fltarr(z_i.length,n_elements(e_spec[0,*]))
    ; set first frequency component to zero (mean of time series)
    z_spec[*,0] = 0
    ; set second component to interpolation
    z_spec[*,1] = z_i
    ; loop to interpolate rest of frequencies
    for w=2, n_elements(z_spec[0,*])-1 do begin
      z_d = alog10(e_spec[*,w])
      z_i = griddata(x_d,y_d,z_d, method='Kriging', xout=x_grid,yout=y_grid,variogram=[1]) 
      z_spec[*,w] = z_i
    endfor
    
    ; create time series of data
    ; from the interopolate z spectrum
    ts_ran = ts_random(reform(10.^z_spec[0,*]), d_t, r_s=long(r_seed))
    psd = psd_calc(ts_ran, res=d_t)
    
    ts_arr = fltarr(z_i.length,ts_ran.length)
    ts_norm = ts_arr
    
    ; get a set of random numbers to generate the time series
    ; increment the seed by i each time so each hour we generate
    ; a new set of seeds different from the last hour
    r_ts = randomn(long(r_seed)+i, z_i.length)
    r_ts = r_ts-min(r_ts)
    r_ts = long(r_ts*1000.)
    for w=0L, z_i.length-1 do begin 
      ts = ts_random(reform(10.^z_spec[w,*]), d_t, r_s=long(r_ts[w]))
      
      ts_arr[w,*] = ts
      ts_norm[w,*] = normalize_vec(ts)
    endfor  
    
    ;add the new timeseries to the final list
    tf.add, ts_arr 
    
    ; sum data and interpolated arrays
    ; over all frequencies for visulaization
    dat_p = total(alog10(e_spec),2,/nan)
    int_p = total(z_spec,2,/nan)
  
    c_min = min([dat_p,int_p],max=c_max)
  
    if verbose gt 0 then begin
    
      fixplot
      window,1, xsize = 500, ysize = 1000
      !p.multi=[0,2,6,0,0]
      !p.charsize=1.5
      !y.omargin = [5,5]
      !y.margin = [0,0]
      
      ;plot data
      plot, x_grid, y_grid, /isotropic, /nodata, title='Mapped Summed PSD!DE!N', $
        xtickf='neg_ticks', ytickf='neg_ticks'
      loadct,25,/silent
      ;cc = bytscl(alog10(total(e_spec,2))
      plots, x_d, y_d, color=bytscl(dat_p,min=cmin,max=cmax), psym=sym(1) 
      
      ; plot interpolated data
      loadct,0,/silent
      plot, x_grid, y_grid, /isotropic, /nodata, title='Interpolated Summed PSD!DE!N', $
        xtickf='neg_ticks', ytickf='neg_ticks'
      loadct,25,/silent
      plots, x_grid, y_grid, color=bytscl(int_p,min=cmin,max=cmax), psym=sym(1)
      
      loadct,0,/silent
      xyouts,0.5,0.98, time_string(t_arr[i])+'-'+time_string(t_arr[i]+3600,tformat='hh:mm:ss'), $
         /normal, alignment=0.5, charsize=1.
      
      ;get a random set of 
      !p.multi=[5,1,6,0,0]
      !p.charsize=1.5
      !y.margin = [0,5]
      !x.margin = [10,10]
      t = findgen(n_elements(ts_arr[0,*]))*d_t
      t_range = [min(t,max=max_t),max_t]
      f_range = [0,max(e_map.freq)]
      y_range = [0,z_i.length-1]
      
      loadct,0,/silent
      plot, o1h.date_themis, o1h.dst , _extra=x_tk, ytitle='Dst'
      oplot, [t_arr[i],t_arr[i]],!y.crange,linestyle=2
      oplot, [t_arr[i],t_arr[i]]+3600,!y.crange,linestyle=2 
      
      ; plot interpolated spectrum
      loadct,0,/silent
      plot, f_range,y_range,/nodata,ytitle='Interpolated PSD'
      loadct,25,/silent
      tvscale,transpose(z_spec),/overplot,/nointerpolation
      
      ; plot random time series
      loadct,0,/silent
      plot, t,y_range,/nodata, ytitle='Normalized!CRandom TS', title='Random TS, '+strtrim(long(d_t),2)+' sec resolution'
      loadct,25,/silent
      tvscale,transpose(ts_norm),/overplot,/nointerpolation
  
      !p.multi=[1,1,3,0,0]
      !p.charsize=1.5
      !y.margin = [0,5]
      !x.margin = [15,20]
      ts_num = normalize_vec((randomn(ssss,1000)+1))
      ts_num = abs(ts_num[0:10])*(z_i.length-1)
      ts_num = long(ts_num)
      l_leg = sqrt(x_grid[ts_num]^2. + y_grid[ts_num]^2.)
      m_leg = atan(y_grid[ts_num],x_grid[ts_num])*(180./!pi)*(24./360.)
      bd = where(m_leg lt 0, c)
      if c gt 0 then m_leg[bd] = m_leg[bd]+24.
      y_leg = string(l_leg, format='(F4.1)')+', '+string(m_leg,format='(F4.1)')  
      stack_plot,t,transpose(ts_arr[ts_num,*]),y_leg, $
        title='Random Sample of Time Series', xtitle = 'Time - seconds', $
        ytitle='Electric Field Perturbation - mV/m'
      xyouts,!p.clip[2], !p.clip[3], ' L/R, MLT',/device
    
      if keyword_set(movie) then begin
        wset,1
        makepng, o_dir+'ran_ts'
        im_mov = Read_PNG(o_dir+'ran_ts.png')
        void = video.Put(stream, im_mov)
      endif
    endif
    
    
  endfor
  
  if keyword_set(movie) then video.Cleanup

  
  r_grid = sqrt(x_grid^2+y_grid^2.)
  mlt_grid = atan(y_grid,x_grid)*180/!pi * 24/360.
  neg_mlt = where(mlt_grid lt 0)
  mlt_grid[neg_mlt] = mlt_grid[neg_mlt]+24

  ts_final = tf.ToArray() ;in mV/m
  ts_final = ts_final/1000. ;convert to V/m

  if w_dat eq 1 then begin
    h5_file = time_string(sd,tformat='YYYYMMDD')+'_'+time_string(ed,tformat='YYYYMMDD')
    h5_file = 'ULF_ran_ts_'+h5_file
    h5_file = o_dir+h5_file+'.h5'
    
    if (file_test(h5_file)) then file_delete, h5_file
    
    mg_h5_putdata, h5_file, 'tstamp',t_arr
    mg_h5_putdata, h5_file, 'tstamp.attributes','Time of hourly time series, number of seconds since 1970  (UNIX time), units [s]'
    
    mg_h5_putdata, h5_file, 'delta_t',d_t
    mg_h5_putdata, h5_file, 'delta_t.attributes','Time step of each hourly time series
    
    
    mg_h5_putdata, h5_file, 'r_grid',r_grid
    mg_h5_putdata, h5_file, 'r_grid.attribute','Radial grid position, units [RE]'
    
    mg_h5_putdata, h5_file, 'mlt_grid',mlt_grid
    mg_h5_putdata, h5_file, 'mlt_grid.attributes','Magnetic Local Time grid position, units [Hr]'
    
    mg_h5_putdata, h5_file, 'x_grid', x_grid
    mg_h5_putdata, h5_file, 'x_grid.attributes', 'X axis grid position, units [RE]'
     
    mg_h5_putdata, h5_file, 'y_grid', y_grid
    mg_h5_putdata, h5_file, 'y_grid.attributes', 'Y axis grid position, units [RE]'
    
    mg_h5_putdata, h5_file, 'ts_random', ts_final
    mg_h5_putdata, h5_file, 'ts_random.attributes', 'Array of randomized phase time series for the azimuthal electric field, ' + $
      'array is [time,grid,hourly time series], ' + $ 
      'array has dimensions [# elements tstamp, # elements grid, 3600/delta_t], resolution - delta t, units [V/m]'
      
    mg_h5_dump, h5_file
  endif

end



;MAIN
;



;ulf_random,'2016-09-26','2016-10-07',verbose=1, /movie

ulf_random,'2016-09-26','2016-09-27',verbose=0, w_dat=1

end