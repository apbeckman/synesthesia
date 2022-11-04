vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// Aya Tunnel by Martijn Steinrucken aka BigWings - 2019
// countfrolic@gmail.com
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// 
//
// Inspired by a vision I had during an ayahuasca ceremony.
// This is a 'sweetend up' version of the tunnel I felt I was flying through.
//
// This is just 3 octaves of gyroids
// First two layers are intersected, next two are mixed.
//

// Music: Asura - Crossroads Limiter
// https://soundcloud.com/licoknigi-1/asura-crossroads-limiter


#define SURF_DIST .001
#define MAX_DIST 75.
#define MAX_STEPS 300

// use these to tweak the effect
#define KALEIDOSCOPE
#define PERLESCENCE
#define SHADOW

vec3 lightPos;

float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float GetHeight(float depth) {
	return 1.5-sin(depth*.5);
}

vec3 GetRayDir(vec2 uv, vec3 ro, vec3 lookat, vec3 up, float zoom) {
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(up, f)),
        u = cross(f, r),
        c = ro + f*zoom,
        i = c + uv.x*r + uv.y*u,
        rd = normalize(i-ro);
    return rd;
}

float LineDist(vec3 p, vec3 a, vec3 b) {
	vec3 ba = a-b,
        pa = a-p;
    
    return length(cross(ba, pa))/length(ba);
}

mat2 Rot(float a) {
    float s = sin(a), c = cos(a);
    return mat2(c, -s, s, c);
}

float GetDist(vec3 p) {
	vec2 m = iMouse.xy/RENDERSIZE.xy - .5;

   float lightDist = length(p-lightPos);
    
    p.z += smoothTime*0.75;
    float z = p.z*.5;
     
    p.y += sin(z);
    
    
    float a = smoothstep(0.,3., lightDist);
    float w = mix(-1., 1., sqrt(a));
    float d = (abs(dot(sin(p), cos(p.zxy)))-.1*w);
    
    float s = 5.;
    p.x += smoothTime*.2;
    
    float second = (abs(dot(sin(p*s), cos(p.zxy*s)))-w)/s;
    
    s = 8.;
    float second2 = (abs(dot(sin(p*s), cos(p.zxy*s)))-w)/s;
    
    second = mix(second, second2, sin(z*.2)*.5+.5);
    
    d = smin(d, second,-.1);
    
    
    m.y=sin(p.z*.345)*.5-.2;
    s = 26.2;
    float third = (abs(dot(sin(p*s), cos(p.zxy*s)))-.1)/s;
    d = mix(d, third, m.y);
    
    return d*.25;
}

vec3 GetNormal(vec3 p) {
    vec2 e = vec2(.01, 0);
    vec3 n = GetDist(p) - vec3(
        GetDist(p-e.xyy),
        GetDist(p-e.yxy),
        GetDist(p-e.yyx)
    );
    
    return normalize(n);
}

vec4 CastRay(vec3 ro, vec3 rd) {
    float dS, dO;
    vec3 p;
    for(int i=0; i<MAX_STEPS; i++) {
    	p = ro + dO * rd;
        dS = GetDist(p);
        dO += dS;
        
        if(dS<SURF_DIST || dO>MAX_DIST) break;
    }
    
    return vec4(p, dS);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec2 UV = uv;
    vec3 col = vec3(0);
	
    float t = smoothTime*0.5;
    float pulse = sin(smoothTimeC*.125)<-0.5?1.:0.;
    vec2 m = iMouse.xy/RENDERSIZE.xy;
    
    #ifdef KALEIDOSCOPE
    uv /= dot(uv,uv)*2.;
    uv = abs(uv);
    float kaleidoscope = smoothstep(.4, .6, sin(-t*.2));
    uv = mix(UV, uv, kaleidoscope);
    #endif
    
    float z = t*.5;
    
    float zoom = 1.;
    vec3 up = vec3(sin(t*.2),cos(t*.2),0);
    vec3 ro = vec3(0, GetHeight(t), 0);
    vec3 lookat = vec3(0,GetHeight(t+.5),1);
    vec3 rd = GetRayDir(uv, ro, lookat, up, zoom);
    float lightDepth = 7.+sin(t*.162-2.)*-2.;
    lightPos = vec3(0,GetHeight(t+lightDepth),lightDepth);

    
    vec4 hit = CastRay(ro, rd);
    vec3 p = hit.xyz;
    float dS = hit.w;
    float dO = length(ro-p);
    
    
    vec3 glowCol = sin(vec3(.123, .234, .345)*smoothTimeB*0.25)*.5+.85;
    col = glowCol;
    
    vec3 L = lightPos-p;	// vec from surface to light
    float ld = length(L);	// distance from surface to light


    vec3 n = GetNormal(p);

    vec3 l = L/ld;

    float dif = clamp(dot(n,l), 0., 1.);
    float att = 1./ld;

    float trust = smoothstep(.8, 1., sin(t))*.6+.4;
    //trust = 1.;

    col = dif*att*trust*glowCol*40.;
    col *= dif*att*trust;

    float spec = clamp(dot(reflect(rd, n), l),0.,1.);
    spec = pow(spec, 20.);
    col += spec*trust*att;

    /*col=1.-col;
col *= .23;
col *= vec3(.8,.7,.9);*/

    if(length(lightPos-ro)<dO) {
        float minLd = LineDist(lightPos, ro, ro+rd); 
        vec3 light = glowCol*trust*.25/minLd;
        col += light;
    }
    col *= smoothstep(5., 6.5, t);

    #ifdef SHADOW
    vec4 shadow = CastRay(lightPos, -l);
    if(abs(length(shadow.xyz-lightPos)-ld) > SURF_DIST) {
        col *= .05;
    }
    #endif

    float fogPhase = sin(t*.1);
    fogPhase *= fogPhase*fogPhase;

    col = mix(col, glowCol, fogPhase*length(hit.xyz-ro)/20.);

    #ifdef PERLESCENCE
    vec3 perl = n*.5+.5;
    perl*=perl;
    pulse = pow(sin(-ld*.4+t)*.5+.5,106.);

    perl = sin(perl*t*.005)*.5+.5;
    vec3 perlCol = pulse*perl*mix(glowCol,vec3(1), 2.75);

    pulse = sin(t*-.12)*.5+.5;
    pulse = pulse*pulse;
    //col = mix(col, perlCol, pulse);

    col += perlCol*pulse;
    #endif


   
  // if(UV.x>m.x*1.77-.5)col = 1.-col;
   
    #ifdef KALEIDOSCOPE
    col *= 1.+kaleidoscope;
    
    #endif
    col = pow(col, vec3(sin(t*.0123)*.25+1.5));
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}