pro niceprint, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14

if (keyword_set(v14)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i),v12(i),v13(i),v14(i)
    return
endif 
if (keyword_set(v13)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i),v12(i),v13(i)
    return
endif 
if (keyword_set(v12)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i),v12(i)
    return
endif 
if (keyword_set(v11)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i),v10(i),v11(i)
    return
endif 
if (keyword_set(v10)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i),v10(i)
    return
endif 
if (keyword_set(v9)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i),v9(i)
    return
endif 
if (keyword_set(v8)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i),v8(i)
    return
endif 
if (keyword_set(v7)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i),v7(i)
    return
endif
if (keyword_set(v6)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i),v6(i)
    return
endif
if (keyword_set(v5)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i),v5(i)
    return
endif
if (keyword_set(v4)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i),v4(i)
    return
endif
if (keyword_set(v3)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i),v3(i)
    return
endif
if (keyword_set(v2)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i),v2(i)
    return
endif
if (keyword_set(v1)) then begin 
    for i=0,n_elements(v1)-1 do $
      print,v1(i)
    return
endif

return
end
