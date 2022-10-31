

			//******** BuffA Code Begins ********

// Thanks to dharrys for the noise reduction code!
// Source: https://www.shadertoy.com/view/MtsBWj
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord / RENDERSIZE.xy;
	float offset = 1.0 / RENDERSIZE.x;
	vec3 o = vec3(-offset, 0.0, offset);
	vec4 gx = texture(syn_UserImage, uv + o.xz);
	vec4 gy = gx;
	vec4 t;

    //Calculate sobel
	gx += 2.0*texture(syn_UserImage, uv + o.xy);
	t = texture(syn_UserImage, uv + o.xx);
	gx += t;
	gy -= t;
	gy += 2.0*texture(syn_UserImage, uv + o.yz);
	gy -= 2.0*texture(syn_UserImage, uv + o.yx);
	t = texture(syn_UserImage, uv + o.zz);
	gx -= t;
	gy += t;
	gx -= 2.0*texture(syn_UserImage, uv + o.zy);
	t = texture(syn_UserImage, uv + o.zx);
	gx -= t;
	gy -= t;
    
	vec4 sobel = sqrt(gx*gx + gy*gy);
    vec4 blur_=texture( syn_UserImage, vec2(uv.x,uv.y), 1.5);
    
    fragColor = mix(blur_, texture(syn_UserImage, uv), sobel);
	return fragColor; 
 } 


#define MAX_ITER 256
#define MAX_DIST 100.
#define MIN_DIST .01

float luma(vec3 v) {
    return .2126*v.x + .7152*v.y + .0722*v.z;
}

vec2 wallUV(vec3 p) {
    float aspect = RENDERSIZE.x/RENDERSIZE.y;
    return (p.xy*vec2(1./aspect,1.))/7.5+vec2(.5,.5);
}

float dstPlane(vec3 p, vec4 n) {
    return dot(p,n.xyz) + n.w;
}

float dstScene(vec3 p) {
    float disp = luma(texture(BuffA, wallUV(p)).xyz);
    return dstPlane(p, vec4(0.,0.,-1.,2.)) + disp * .5;
}

float raymarch(vec3 ro, vec3 rd) {
    float t = 0.;
    for(int i = 0; i < MAX_ITER; i++) {
        float d = dstScene(ro+rd*t);
        if(d < MIN_DIST || d > MAX_DIST) {
            break;
        }
        t += d * .75;
    }
    return t;
}

// Shadows by iq
// Source: https://www.shadertoy.com/view/Xds3zN
float softshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax, in float hardness )
{
	float res = 1.0;
    float t = mint;
    for( int i=0; i<32; i++ )
    {
		float h = dstScene( ro + rd*t );
        res = min( res, hardness*h/t );
        t += clamp( h, 0.06, 0.30 );
        if( h<0.001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

vec3 calcNormal(vec3 p, float t) {
    vec2 e = vec2(t * MIN_DIST, 0.);
    vec3 n = vec3(dstScene(p+e.xyy)-dstScene(p-e.xyy),
                  dstScene(p+e.yxy)-dstScene(p-e.yxy),
                  dstScene(p+e.yyx)-dstScene(p-e.yyx));
    return normalize(n);
}

vec3 pixel(vec3 ro, vec3 rd) {
    float t = raymarch(ro, rd);
    if(t < MAX_DIST) {
        vec3 p  = ro + rd * t;
        vec3 n  = calcNormal(p, t);
        vec3 r  = normalize(reflect(rd, n));
        vec2 uv = wallUV(p);
        
        vec3 mat = texture(image46, uv).xyz;
        vec3  l = normalize(vec3(45.,30.,-45.));
        float a = softshadow(p+l*MIN_DIST,l,MIN_DIST+.05,MAX_DIST, 3.);
        float d = max(dot(l,n)*a,.2);
        float s = pow(max(dot(l,r),0.),5.+25.*luma(mat))*a;
        
        return mat*d+s;
    }
    return rd;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord-RENDERSIZE.xy*.5)/RENDERSIZE.y;
    vec3 ro = vec3(0., 0., -5.);
    vec3 rd = vec3(uv, 1.);
    
    fragColor.xyz = pixel(ro, normalize(rd));
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}