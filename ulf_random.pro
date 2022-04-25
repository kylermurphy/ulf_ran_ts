pro ulf_random, sdate, edate, res, $
  l_min=l_min, $
  l_max=l_max, $
  l_bin=l_bin, $
  mlt_min=mlt_min, $
  mlt_max=mlt_max, $
  mlt_bin=mlt_bin, $
  d_dir = d_dir
  
  if keyword_set(d_dir) then d_dir=d_dir else d_dir='D:\data\magnetometer\psd\'
  
  if res eq 1 then begin
    npts = 1801
  endif else if res eq 10 then begin
    npts = 181
  endif else begin
    message, 'Resolution must be 10 or 1'
  endelse
  
  ; create an hourly time series 
  ; for readin in the PSD files
  sd = time_double(sdate)
  ed = time_double(edate)
  t_arr = sd+dindgen((ed-sd)/3600.)*3600.+1
  t_arr = time_string(t_arr, tformat='YYYY-MM-DD/hh')
  t_arr = time_double(t_arr)
  
  y_i = median(long(time_string(t_arr, tformat='YYYY')))
  stn_vals = stn_dat(y_i)
  
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
    spec = fltarr(stn.length, npts)
    
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
      
      spec[r_c,*] = r_s.psd[h_i,0:npts-1]
      lsh[r_c] = stn_vals.lshell[s_p]
      mlt[r_c] = m_val
      
      print, strtrim(h_i,2)+' '+stn[j]+' '+strtrim(m_val,2)+' '+strtrim(lsh[r_c])
      r_c++
    endfor
  
    spec = spec[0:r_c-1,*]
    lsh = lsh[0:r_c-1]
    mlt = mlt[0:r_c-1]
  
    stop
  endfor
  
  stop




end



;MAIN
;



ulf_random,'2016-09-26','2016-10-07',10

end