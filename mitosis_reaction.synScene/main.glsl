
int media_switch = int(_isMediaActive());
vec4 media = mix(vec4(0.0), ((texture(media_pass, _uv))), media_switch);

vec4 mediaEdges = mix(vec4(0.0), texture(media_pass_fx, _uv), media_switch);
vec4 edgeMix = mix(media, mediaEdges, 0.125);
vec4 media_mask = edgeMix + _invert(_loadMediaAsMask(), 1.0);
float media_lum = PI *.25 * sin(length(edgeMix));
float inv_med_lum = sin(length(_invert(edgeMix, 1.0)));
#define R RENDERSIZE.xy
float a = 2.0;
// #define a 2.
float o(float a){return a * .7071;} 

// #define o2 a2*.7071
vec4 A (vec2 U) {
    return texture(BuffA,U/R);
}

vec4 B (vec2 U) {
    return texture(BuffB,U/R);
}

vec4 C (vec2 U) {
    return texture(BuffC,U/R);
}

vec4 D (vec2 U) {
    return texture(BuffD,U/R);
}

float smin(float c, float b, float k) {
    float h = max(k - abs(c - b), 0.) / k;
    return min(c, b) - h * h * h * k * 1. / 6.;
}

float ln (vec3 p, vec3 aa, vec3 b) {return length(p-aa-(b-aa)*dot(p-aa,b-aa)/dot(b-aa,b-aa));}


float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d = dist(RENDERSIZE.xy*0.5,_xy.xy)*(_mouse.x/RENDERSIZE.x+0.1)*0.005;

//saving some space and making mods easier by creating a single method to call after the render pass-specific buffer operations
vec4 fx_body(vec4 Q, vec4 n) {
    vec2 U = _xy;
    n *= .125;
    
    vec4 dx = n-Q;

    Q += dx*vec4(1.,.3,1.,1);
    float x = .3*Q.x*Q.y*(1.-Q.y); 
    x += growth*0.0125;
    x += pow(syn_BassLevel*0.75, 2.0)*0.001;
    x -= media_lum *0.005;
    // Q.y = Q.y+x-0.025-.1*Q.z;
    Q.y = Q.y+x-0.025-.1*Q.z;
    Q.y += growth *0.01;
    // Q.x += media_lum*0.075;
    // Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    Q.y += media_lum * 0.0125;
    Q.z = Q.z*0.98+ .6*dx.y+(.001)*Q.y + 0.01* split;
    Q.z += media_lum * 0.00025;

    Q = clamp(Q,0.,1.);

    if (_mouse.z> 0. && length(U-_mouse.xy) < smin(30., 40., 0.2)) Q.y=1.;
    
    if (FRAMECOUNT <= 1) {
        Q = vec4(1,0,0,0);
        if (length(U-0.5*R)<10.) Q.y += 1.;
        if (length(U-0.45*R)<10.) Q.y += 1.;

    }
	return Q; 
}

			//******** BuffA Code Begins ********


vec4 renderPassA() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    Q = D(U);
    float o = o(a);
    float aa = TIME+3.*length(U), c = cos(aa), s = sin(aa);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = D(U+vec2(0,a)*m)+D(U+vec2(a,0)*m)+D(U+vec2(0,-a)*m)+D(U+vec2(-a,0)*m)+D(U+vec2(-o,o)*m)+D(U+vec2(o,-o)*m)+D(U+vec2(-o,-o)*m)+D(U+vec2(o,o)*m);
    // n *= .125;
    // vec4 dx = n-Q;
    // Q += dx*vec4(1.,.3,1.,1);
    // float x = .3*Q.x*Q.y*(1.-Q.y);
    // Q.y = Q.y+x-0.025-.1*Q.z;
    // Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    // Q = clamp(Q,0.,1.);
    // if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    // if (FRAMECOUNT < 1) {
    //     Q = vec4(1,0,0,0);
    //     if (length(U-0.5*R)<10.) Q.y += 1.;
    // }
    Q = fx_body(Q, n);
    Q.y += media_lum*0.0125;
    Q.y += media_lum*0.0125;
	return Q; 
 } 


			//******** BuffB Code Begins ********

#define R RENDERSIZE.xy
vec4 renderPassB() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    float a1 = a + 1.0;
    float o = o(a);
    Q = A(U);
    float aa = TIME+3.*length(U), c = cos(aa), s = sin(aa);
    mat2 m = mat2(c,-s,s,c+sin(TIME));
    vec4 n = A(U+vec2(0,a1)*m)+A(U+vec2(a1,0)*m)+A(U+vec2(0,-a1)*m)+A(U+vec2(-a1,0)*m)+A(U+vec2(-o,o)*m)+A(U+vec2(o,-o)*m)+A(U+vec2(-o,-o)*m)+A(U+vec2(o,o)*m);
    // n *= .125;
    // vec4 dx = n-Q;
    // Q += dx*vec4(1.,.3,1.,1);
    // float x = .3*Q.x*Q.y*(1.-Q.y);
    // Q.y = Q.y+x-0.025-.1*Q.z;
    // Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    // Q = clamp(Q,0.,1.);
    // if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    // if (FRAMECOUNT <= 1) {
    //     Q = vec4(1,0,0,0);
    //     if (length(U-0.5*R)<10.) Q.y += 1.;
    // }
    Q = fx_body(Q, n);
    return Q; 
 } 


			//******** BuffC Code Begins ********

#define R RENDERSIZE.xy
vec4 renderPassC() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    float a1 = a + 1.0;
    float o = o(a);
    Q = B(U);
    float aa = TIME+3.*length(U), c = cos(aa), s = sin(aa);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = B(U+vec2(0,a1)*m)+B(U+vec2(a1,0)*m)+B(U+vec2(0,-a1)*m)+B(U+vec2(-a1,0)*m)+B(U+vec2(-o,o)*m)+B(U+vec2(o,-o)*m)+B(U+vec2(-o,-o)*m)+B(U+vec2(o,o)*m);
    // n *= .125;
    // vec4 dx = n-Q;
    // Q += dx*vec4(1.,.3,1.,1);
    // float x = .3*Q.x*Q.y*(1.-Q.y);
    // Q.y = Q.y+x-0.025-.1*Q.z;
    // Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    // Q = clamp(Q,0.,1.);
    // if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    // if (FRAMECOUNT < 1) {
    //     Q = vec4(1,0,0,0);
    //     if (length(U-0.5*R)<10.) Q.y += 1.;
    // }
    Q = fx_body(Q, n);
    Q.y += media_lum*0.005;
    Q.y += media_lum*0.005;
	return Q; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;
    float a2 = a + 2.0;
    float o = o(a);
    Q = C(U);
    float aa = TIME+3.*length(U), c = cos(aa), s = sin(aa);
    mat2 m = mat2(c,-s,s,c);
    vec4 n = C(U+vec2(0,a2)*m)+C(U+vec2(a2,0)*m)+C(U+vec2(0,-a2)*m)+C(U+vec2(-a2,0)*m)+C(U+vec2(-o,o)*m)+C(U+vec2(o,-o)*m)+C(U+vec2(-o,-o)*m)+C(U+vec2(o,o)*m);
    // n *= .125;
    // vec4 dx = n-Q;
    // Q += dx*vec4(1.,.3,1.,1);
    // float x = .3*Q.x*Q.y*(1.-Q.y);
    // Q.y = Q.y+x-0.025-.1*Q.z;
    // Q.x = Q.x+0.1*Q.x*(1.-Q.x)-x;
    // Q.z = Q.z*0.98+ .6*dx.y+.001*Q.y;
    // Q = clamp(Q,0.,1.);
    // if (_mouse.z> 0. && length(U-_mouse.xy) < 40.) Q.y=1.;
    // if (FRAMECOUNT <= 1) {
    //     Q = vec4(1,0,0,0);
    //     if (length(U-0.55*R)<10.) Q.y += 1.;
    //     if (length(U-0.45*R)<10.) Q.y += 1.;
    // }
    Q = fx_body(Q, n);
    Q.y += media_lum*0.00125;
    Q.y += media_lum*0.00125;
	return Q; 
 } 



vec4 renderMainImage() {
	vec4 Q = vec4(0.0);
	vec2 U = _xy;

    vec4
        n = A(U+vec2(0,1)),
        e = A(U+vec2(1,0)),
        s = A(U-vec2(0,1)),
        w = A(U-vec2(1,0));

    vec2 h = vec2(e.w-w.w,n.w-s.w);
    vec2 g = vec2(e.y-w.y,n.y-s.y);
    vec4 b = A(U);
   	vec4 dx = 0.25*(n+e+s+w) - b;

	Q = (b.w+b.y)*normalize(abs(sin(4.+.5*(b.z-b.x)*vec4(3,2,1,4))));
    
    vec3 no = normalize(vec3(g-h,-1));
    vec3 r = reflect(vec3(0,0,1),no);
    vec3 l = vec3(2.*R.x,2.*R.y,10.*R.y);
    vec3 u = vec3(U,0);
    vec3 lu = l-u;
    float oo = length(r-lu*dot(r,lu)/dot(lu,lu));
    Q = (Q)*(exp(-oo)+2.*exp(-7.5*oo)+1e2*exp(-30.*oo));
    // Q *= 1.+texture(iChannel1,r);
    no = normalize(vec3(g-h,1));
    r = refract(vec3(0,0,1),no,.2);
    l = vec3(3.*R.x,3.*R.y,10.*R.y);
    u = vec3(U,0);
    lu = l-u;
    oo = length(r-lu*dot(r,lu)/dot(lu,lu));
    Q = 2.*dx+(Q*0.9+0.1)*(0.4+exp(-oo)+10.*exp(-9.5*oo));
    // Q *= 1.+0.2*sin(Q+3.*texture(iChannel1,r).x*vec4(1,2,3,4));
    // Q *= 1.0 + media_lum*0.1;
    Q += mix(vec4(0.0), media_mask, see_media) * 0.1;

	return Q; 
 } 

vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _loadMedia();
    return media;
}

vec4 mediaPassFX() {
    vec4 media_fx = _edgeDetectSobel((media_pass));
    return media_fx;
}

vec4 renderMain() {
    if(PASSINDEX == 0) {
        return renderPassA();
    }
    if(PASSINDEX == 1) {
        return renderPassB();
    }
    if(PASSINDEX == 2) {
        return renderPassC();
    }
    if(PASSINDEX == 3) {
        return renderPassD();
    }
    if(PASSINDEX == 4) {
        return mediaPass();
    }
    if(PASSINDEX == 5) {
        return mediaPassFX();
    }
    if(PASSINDEX == 6) {
        return renderMainImage();
    }

}