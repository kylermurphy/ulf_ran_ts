function ts_random, psd, res, r_s=r_s

  if keyword_set(r_s) then r_s = r_s else r_s = 17

  if (psd.length ne 1801 and psd.length ne 181) then message,'Spectral length must be 1801 or 181 and derived from 1 hour time series'
  
  ts_len = psd.length*2-2
  
  ; convert psd to 1-sided FFT spectrum
  f_spec = sqrt(psd/(res*ts_len*2.))  

  ; create 2-sided FFT from the 1-side FFT
  rfft = fltarr(f_spec.length*2-2)
  rfft[0:f_spec.length-1] = f_spec
  rfft[f_spec.length,-1] = reverse(f_spec[1:-2])

  ; create random phase vector
  rphi = fltarr(f_spec.length*2-2)
  rphi[1:ts_len/2] = randomn(long(r_s), n_elements(rphi[1:ts_len/2]))*2*!pi
  rphi[ts_len/2+1:-1] = -1*reverse(rphi[1:ts_len/2-1])
  
  ; combine phase and fft to generate fft
  rp_fft = rfft*exp(complex(0,1)*rphi)
  ; calculate the inverse fft to 
  ; get a time series back
  rp_ifft = fft(rp_fft,/inverse)
  ; take the real part of the inverse
  ; fft 
  ; complex part should be close to 
  ; zero
  r_ts = real_part(rp_ifft)
  
  return, r_ts
  
end
