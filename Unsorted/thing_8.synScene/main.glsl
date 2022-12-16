

			//******** Common Code Begins ********

vec2 R;
vec4 M;
int I;
#define N 9.
#define R3D vec3(R/N,N*N)
#define e(a) mat2(cos(a),-sin(a),sin(a),cos(a))
#define d2(U) ((U).xy+vec2(mod(floor((U).z),N),floor(floor((U).z)/N))*R/N)
#define d3(u) vec3(mod(u,R/N),floor(u/R*N).x+floor(u/R*N).y*N)
#define _3D  vec3 U = d3(u)
#define Sampler vec4 T(vec3 U) {return mix(texture(iChannel0,d2(vec3(U.xy,floor(U.z)))/R),texture(iChannel0,d2(vec3(U.xy, ceil(U.z)))/R),fract(U.z));}
#define A(U) texture(cha,d2(mod(U,R3D))/R)
#define B(U) texture(chb,d2(mod(U,R3D))/R)
#define Main void mainImage (out vec4 Q, in vec2 U)
float signe (float x) {return atan(100.*x);}
void prog (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = vec4(0); b = vec4(0);
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
    	vec4 aa = A(U+u), bb = B(U+u);
        aa.xyz += bb.xyz;
        #define q 1.1
		vec3 w1 = clamp(aa.xyz-0.5*q,U - 0.5,U + 0.5),
             w2 = clamp(aa.xyz+0.5*q,U - 0.5,U + 0.5);
        float m = (w2.x-w1.x)*(w2.y-w1.y)*(w2.z-w1.z)/(q*q*q);
        aa.xyz = 0.5*(w1+w2);
        a.xyz += aa.xyz*aa.w*m;
        b.xyz += bb.xyz*aa.w*m;
        a.w += aa.w*m;
    }
    if (a.w>0.) {
        a.xyz/=a.w;
        b.xyz/=a.w;
    }
}
void prog2 (vec3 U, out vec4 a, out vec4 b, sampler2D cha, sampler2D chb) {
	
    a = A(U); b = B(U);
    vec3 f = vec3(0); float m = 0.;
    for (int x = -1; x <= 1; x++)
    for (int y = -1; y <= 1; y++)
    for (int z = -1; z <= 1; z++)
    {
        vec3 u = vec3(x,y,z);
        float l = length(u);
        if (l>0.) {
    		vec4 aa = A(U+u), bb = B(U+u);
            f += 1e-2*(aa.w*(1.-.2*aa.w))*u/l;
            m += aa.w;
        }
    }
    if (m>0.) b.xyz += f/m;
    
    
    // Boundaries:
   	b.xyz -= 1e-3*signe(a.w)*(a.xyz-0.5*R3D)*sin(1e-5*float(I));

    
    if (I<1||U.x<1.||R3D.x-U.x<1.||R3D.y-U.y<1.||R3D.x-U.x<1.||U.z<1.||R3D.z-U.z<1.) {
    	a = vec4(U,0);
        b = vec4(0);
        if (length(U-0.5*R3D) < 0.2*R3D.y) a.w = 20.;
    }
}

			//******** BuffA Code Begins ********

vec2	_TexelSize;
float	_Scale;
float	_Tresh;

//#define FEED_DEFAULT_ .0550
//#define KILL_DEFAULT .0620

//#define FEED_MITOSIS .0367
//#define KILL_MITOSIS .0649

vec2	laplacian_convolution(sampler2D tex, vec2 uv, vec2 texsize)
{
	vec2	ret = vec2(0.);
    
    if (uv.x == 0. || uv.y == 0. || uv.x== 1. || uv.y ==1.)
        return (ret);
    ret += texture(tex, vec2(uv.x , uv.y) ).xy * -1.;
    
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x , uv.y -texsize.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x , uv.y +texsize.y) ).xy * (.2);
    
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y -texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y -texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y +texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y +texsize.y) ).xy * (.05);
    return (ret);
}

vec4 gradient(sampler2D tex, vec2 uv, vec2 offset)
{
        return (texture(tex, uv + offset)) - (texture(tex, uv - offset));
}

vec4 optical_flow(vec2 f)
{
    if (f.x == 0.0 || f.y == 0.0
        || f.x >= RENDERSIZE.x || f.y >= RENDERSIZE.y)
        discard;
    _TexelSize = 1./RENDERSIZE.xy;
    vec2 uv =    f /RENDERSIZE.xy;
	vec4 current = texture(syn_UserImage, uv);
    _TexelSize = 1./RENDERSIZE.xy;
    uv = (f - .0*RENDERSIZE.xy)/RENDERSIZE.xy;
	vec4 prev = texture(BuffA, uv);
    _TexelSize = 1./RENDERSIZE.xy;
    vec2 dx = vec2(_TexelSize.x, 0);
    vec2 dy = vec2(0, _TexelSize.y);
    _TexelSize = 1./RENDERSIZE.xy;
    vec2 ddx = vec2(_TexelSize.x, 0);
    vec2 ddy = vec2(0, _TexelSize.y);

    float	key = ( (current.y) < length(current.xz) ) ? 1. : .0;
    current *= key;
    
    vec4 diff = (current - prev);

    vec4 gx = gradient(BuffA, uv, ddx) + gradient(syn_UserImage, uv, dx);
    vec4 gy = gradient(BuffA, uv, ddy) + gradient(syn_UserImage, uv, dy);

    vec4 gmag = sqrt( abs(gx * gx + gy * gy) + 1.000);
    vec4 invGmag = 1.0 / gmag;
    vec4 vx = diff * (gx * invGmag);
    vec4 vy = diff * (gy * invGmag);

    vec2 flow = vec2(0, 0);
    const float inv3 = 0.33333;
    flow.x = (vx.x + vx.y + vx.z) * inv3;
    flow.y = (vy.x + vy.y + vy.z) * inv3;

    float w = length(flow);
    float nw = (w - _Tresh) / (1.0 - _Tresh);
    flow = mix(vec2(0., 0.), normalize(flow) * nw * _Scale , 1.*step(_Tresh, w));
    return vec4(clamp(flow, .0, 1.) , 0., 1.);
}

vec4 renderPassA() {
	vec4 o = vec4(0.0);
	vec2 f = _xy;

    vec2 R = RENDERSIZE.xy;
    _TexelSize = 1./R;
    _Scale = 1.;
    _Tresh = .003075;
    o = optical_flow(f);
    o.z -= 0.016667;
    o.z += o.y;
    o*= step(_Tresh, length(o) )*1.42;
    //o*=.0;
	return o; 
 } 



			//******** BuffB Code Begins ********

#define FEED_DEFAULT .01 // .0025 // .01 //.01003550
#define KILL_DEFAULT .04 // .014 // .03 //.040620

#define FEED_MITOSIS .0367
#define KILL_MITOSIS .0649
/*
vec2	laplacian_convolution(sampler2D tex, vec2 uv, vec2 texsize)
{
	vec2	ret = vec2(0.);
    
    //if (uv.x == 0. || uv.y == 0. || uv.x== 1. || uv.y ==1.)
    //    return (ret);
    ret += texture(tex, vec2(uv.x , uv.y) ).xy * -1.;
    
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x , uv.y -texsize.y) ).xy * (.2);
    ret += texture(tex, vec2(uv.x , uv.y +texsize.y) ).xy * (.2);
    
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y -texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y -texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x +texsize.x, uv.y +texsize.y) ).xy * (.05);
    ret += texture(tex, vec2(uv.x -texsize.x, uv.y +texsize.y) ).xy * (.05);
    return (ret);
}
*/
vec2 do_life(sampler2D tex, vec2 uv, vec2 texsize)
{
	vec2 ret = vec2(.0);
    
    if (uv.x == 0. || uv.y == 0. || uv.x== 1. || uv.y ==1.)
        return (ret);
    ret += texture(tex, vec2(uv.x , uv.y) ).xy * 1.;
    vec2 sum = vec2(.0);
    
    sum += texture(tex, vec2(uv.x -texsize.x, uv.y) ).xy;
    sum += texture(tex, vec2(uv.x +texsize.x, uv.y) ).xy;
    sum += texture(tex, vec2(uv.x , uv.y -texsize.y) ).xy;
    sum += texture(tex, vec2(uv.x , uv.y +texsize.y) ).xy;
    sum += texture(tex, vec2(uv.x -texsize.x, uv.y -texsize.y) ).xy;
    sum += texture(tex, vec2(uv.x +texsize.x, uv.y -texsize.y) ).xy;
    sum += texture(tex, vec2(uv.x +texsize.x, uv.y +texsize.y) ).xy;
    sum += texture(tex, vec2(uv.x -texsize.x, uv.y +texsize.y) ).xy;
    //sum -= laplacian_convolution(tex, uv, texsize).xy*2.;
    if (ret.x > .0) // alive
    {
    	if (sum.x >= 2. && sum.x <=3.)
        {
            ret.x = 1.;
            ret.y = 1.;
        }
        else
        {
            ret.x = .0;
            ret.y = .0;
        }
    }
    else if (ret.x <= .0)
    {
        if (sum.x == 3.)
        {
            ret.x = 1.;
            ret.y = 1.;
        }
        else
        {
            ret.x = .0;
            ret.y = .0;
        }
    }
    return ret;
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord.xy) / RENDERSIZE.xy;
    vec2 texsize = 1./RENDERSIZE.xy;
    vec4	ret = texture(BuffB, uv);
    
    vec2	ab = ret.xy;
    vec2 mouse = _mouse.xy / RENDERSIZE.xy;
    if (  _mouse.z > 0.0 &&
        mouse.x <= abs(uv.x+.005) && mouse.x >= abs(uv.x-.005) &&
        mouse.y <= abs(uv.y+.005) && mouse.y >= abs(uv.y-.005) )
    {
	    fragColor.x += 1.;
    	ab.y+=1.;
	}
    if (TIME < 2.)
    {
		ab.x = 1.;
	    ab.y = 1.0;
    }
    fragColor *= .0;
    for (float i = .0; i < 1.; i++)
    {//ab*=.995;
      ab.x += texture(BuffA, uv).x;
      ab.y += texture(BuffA, uv).y;
     fragColor.x += clamp(ab.x + (1. * (-.0*1.*laplacian_convolution(BuffA, uv, texsize).x+1.*laplacian_convolution(BuffB, (uv.xy*200.)/200., texsize*1.+texsize*i*.1).x) - ab.x * ab.y * ab.y + FEED_DEFAULT * (1. - ab.x) ) ,0.,1.);
     fragColor.y += clamp(ab.y + (.5 * (-.0*1.*laplacian_convolution(BuffA, uv, texsize).y+1.*laplacian_convolution(BuffB, (uv.xy*200.)/200., texsize*1.+texsize*i*.1).y) + ab.x * ab.y * ab.y - (FEED_DEFAULT + KILL_DEFAULT) * ab.y ),0.,1.);
     uv *= 1.;
     texsize *= 1.;
     //fragColor.xy = do_life(BuffB, uv, texsize);
     //fragColor.xy += .0125*do_life(BuffB, floor(uv*1000.)/1000., texsize*1.).xy;
     if (  _mouse.z > 0.0 &&
        mouse.x <= abs(uv.x+.005) && mouse.x >= abs(uv.x-.005) &&
        mouse.y <= abs(uv.y+.005) && mouse.y >= abs(uv.y-.005) )
	    fragColor.x = 1.;
     //fragColor.x = max(fragColor.x, length(texture(BuffA, uv).xy)*1.);
     //fragColor.x -= .00000001;
         fragColor.xy = clamp(fragColor.xy, .0, 1.);
    }
    fragColor/=1.;
	return fragColor; 
 } 


#define SHOW_DIFUSION 0
#define SHOW_OPTICAL_FLOW 1
#define SHOW_RAYMARCHED_SCENE 2
#define I_MAX	350
#define E		0.00001

//#define BALL
#define SHOW SHOW_RAYMARCHED_SCENE

vec3	cam(vec2 u);
float	map(vec3 p);
float	mylength(vec3 p);
float	mylength(vec2 p);
void	rotate(inout vec2 v, float angle);
vec3	h;
vec4	col;

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 R = RENDERSIZE.xy;
    vec2 uv = (fragCoord.xy - R*.5)/R.yy;
	h = vec3(.0,.0,.0);
    vec4	tex = (texture(BuffB, uv.xy ).xyzy*1.);
    
    col = vec4(.0,.0,.0, 1.);
    
    vec3 pos = vec3(.0, .0, -20.);
    vec3 dir = cam(uv);
    vec3 p;
    vec2 dist = vec2(.0, .0);
    for (float i = .0 ; i < 100.; i++)
    {
        p = pos + dist.y * dir;
        dist.x = map(p);
        dist.y += dist.x;
        if (dist.x < E || dist.y > 150.)
            break;
    }
	h*=.251;

  	#if (SHOW == SHOW_DIFUSION)
	col = texture(BuffB, fragCoord/R);
    #elif SHOW == SHOW_OPTICAL_FLOW
    col = texture(BuffA, fragCoord/R);
    #elif SHOW == SHOW_RAYMARCHED_SCENE
	col.xyz = h;
    #endif
	fragColor.xyz = col.xyz;
    return fragColor;
}

float map(vec3 p)
{
    float mind = 1e5;
    
    vec2 tcd;
    p.z += -15.;
    //rotate(p.xz, -.4+sin(TIME*.5)*.25);
    p.xz = _rotate(p.xz, -.4+sin(TIME*.5)*.25);
    p.zy = _rotate(p.zy,  5.4+TIME*.0);
    p.xz = _rotate(p.xz,  1.+TIME*.0+3.4);
    
    //rotate(p.zy, 5.4+TIME*.0);
    //rotate(p.xz, 1.+TIME*.0+3.4);
    
    vec3 np = normalize(p);
    tcd.x = p.x*(9./RENDERSIZE.x )+.5;
    tcd.y = p.y*(9./RENDERSIZE.y )+.5;//tcd = np.yx;
    #ifdef BALL
    rotate(np.xz, 1.5);
    tcd.x = .5+atan(np.z, np.x)/6.28;
	tcd.y = .5-asin(np.y)/3.14;
    #endif
    tcd = mix(tcd, floor(tcd*300.)/(300.), .5+.5*sin(TIME+(p.x+p.y)*.5) );
    float to = .1*( texture(BuffB, tcd).x*1.);
    
    #ifdef BALL
    mind = min(mind, length(p)-10.+to*10.);
    mind = max(mind, -(length(p)-8.5-to*10. ) );
    #else
    mind = max(p.z-.5, -p.z-.5)-to*10.;
    mind = max(mind, abs(p.x)-16.+.0*-10000.0/RENDERSIZE.x);
    mind = max(mind, abs(p.y)-16.+.0*-10000.0/RENDERSIZE.x);
    #endif
    h += 1.*vec3(.45, .5, .4)*1./max(.1, mind*mind*1. + 15.);
    float mint = max(mind, .01);
    h -= 1.*vec3(.05, .40, .50)*1./max(.00001, mint*mint*1000.001 + 20.);
    
    return mind*.421;
}

vec3 cam(vec2 u)
{
    vec3 up = vec3(.0, 1., .0);
    vec3 ri = vec3(1., 0., .0);
    vec3 fw = vec3(.0, 0., 1.);
	return normalize(up*u.x + ri*u.y + fw);
}

float mylength(vec3 p) {return max(max(abs(p.x), abs(p.y)), abs(p.z));}
float mylength(vec2 p) {return max(abs(p.x), abs(p.y) );}
/*
vec2 rotate(inout vec2 v, float angle)
{
	v = vec2(cos(angle)*v.x+sin(angle)*v.y,-sin(angle)*v.x+cos(angle)*v.y);
	return v; 
 } 
*/


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderMainImage();
	}
}