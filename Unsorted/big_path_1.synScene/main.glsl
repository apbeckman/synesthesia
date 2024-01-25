

#define STEP 256
#define EPS .0001


// from various shader by iq

float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

mat2 m = mat2(.85+Waves,.65+Waves,-.65-Waves,.85-Waves);
 
float noise( in vec2 x )
{
	return sin(1.25*x.x)*sin(1.25*x.y);
}

float fbm6( vec2 p )
{
    p += pow(syn_BassLevel*0.75, 2.)*Twitch;
    float f = 0.0;
    f += 0.500000*(0.5+0.5*noise( p )); p = m*p*2.02;
    f += 0.250000*(0.5+0.5*noise( p )); p = m*p*2.03;
    f += 0.125000*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.062500*(0.5+0.5*noise( p )); p = m*p*2.04;
    f += 0.031250*(0.5+0.5*noise( p )); p = m*p*2.01;
    f += 0.015625*(0.5+0.5*noise( p ));
    return f/(0.96875);
}
vec3 _position;
vec3 camPos = _position + vec3(0.,0.,bass_time * 10.);
float camNoise =  fbm6(camPos.xz * (.025)*1.+0.125*smoothTimeC) *1.;


mat2 getRot(float a)
{
    float sa = sin(a), ca = cos(a);
    return mat2(ca,-sa,sa,ca);
}




float sphere(vec3 center, float radius)
{

    return distance(_position,center+vec3(sin(smoothTime*0.1)*0.5, cos(smoothTime*0.1)*0.5, -sin(smoothTime)*0.125)) - radius;
}

float hozPlane(float height)
{
    return distance(_position.y,height);
}
float idk = pow(syn_BassLevel*0.4+syn_MidLevel*0.4+syn_Level*0.2, 2.);
float swingPlane(float height)
{
    vec3 pos = _position + vec3(0.+sin(smoothTime*0.1)*0.5,0.+cos(smoothTime*0.1)*0.5,bass_time*10.);
    float def =  fbm6((_uvc+pos.xz) * (.125)*1.+0.35*smoothTimeC) * 1.;
    
    float way = pow(abs(pos.x) * (17.*(NoiseStretch)) ,2.5) *.0000125;
    def *= way;

    float ch = height + def;
    
    return max(pos.y - ch,0.);
}

float map(vec3 pos)
{
    _position = pos;
    
    float dist;
    dist = swingPlane(0.);
    
    float sminFactor = 5.25;
    dist = smin(dist,sphere(vec3(0.,-15./Distance,80.*Distance),45.-Distance),sminFactor);
    return dist;
}


vec3 getNormal(vec3 pos)
{
    vec3 nor = vec3(0.);
    vec3 vv = vec3(0.,1.,-1.)*.01;
    nor.x = map(pos + vv.zxx) - map(pos + vv.yxx);
    nor.y = map(pos + vv.xzx) - map(pos + vv.xyx);
    nor.z = map(pos + vv.xxz) - map(pos + vv.xxy);
    nor /= 2.;
    return normalize(nor);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
	vec2 uv = (fragCoord.xy-.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 look = vec3(0.);
    vec3 rayOrigin = vec3(vec2(0.,sin(smoothTime*0.1)*3.+3.0+16.*CamHeight), -1. );
    vec3 rayDir = normalize(vec3(_uvc*FOV , 1.));
    rayDir.xz = _rotate(rayDir.xz, _uvc.x*PI*FOV*Perspective-_uvc.x*PI*Perspective);
    float m_x = smin(rayDir.x*-1, rayDir.x, 0.5);
    float m_y = smin(rayDir.y*-1, rayDir.y, 0.5);
    vec2 rd_mirrored = vec2((m_x), (m_y));
    // vec2 rd_mirrored = vec2(mix(m_x, m_x*-1, invert), mix(m_y, m_y*-1, invert));
    rayDir.xy = mix(mix(rd_mirrored, rd_mirrored*-1, invert) , rayDir.xy, 1.0-vec2(mirror_x, mirror_y));
    rayDir.xz = _rotate(rayDir.xz, PI*Swivel); 
   	rayDir.zy = normalize(getRot(.05*Distance) * rayDir.zy);
   	rayDir.xy = (getRot(.075+-PI*2+sin(bass_time*0.1)*0.0125) * rayDir.xy);
    //rayDir.xy+=_uvc*FOV;
    rayDir.xy = _rotate(rayDir.xy, Spin*PI);
    
    vec3 position = rayOrigin;
    
    
    float curDist;
    int nbStep = 0;
    
    for(; nbStep < STEP;++nbStep)
    {
        curDist = map(position);
        
        if(curDist < EPS)
            break;
        position += rayDir * curDist * .5;
    }
    
    float f;
    //float sound = texture(iChannel0,vec2(.25,.5)).r;
    
    //sound = sin(TIME) * .5 + .5;
    
    float dist = distance(rayOrigin,position);
    f = dist /(99.);
    f = float(nbStep) / float(STEP);
    //removing sound variable and replacing w/my own
    //f += pow(f,3.) * (sound * 75.);
    f *= pow(1.0+highhits*Flash, 2.);
    f *= .8;
    vec3 col = vec3(f);
    
    
    //float shouldColor = 1.- step(f,threshold);
    //col = mix(col,vec3(1.,0.,0.) ,shouldColor);
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}