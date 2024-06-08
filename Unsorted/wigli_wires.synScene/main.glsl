

			//******** BuffA Code Begins ********
#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define C(u) texture(BuffC,(u)/RENDERSIZE.xy)
#define D(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define T(u) texture(image46,(u)/RENDERSIZE.xy)
float lum(vec4 rgb)
{
  return dot(rgb, vec4(0.4, 0.159, 0.11, 0.));
}



/*
#define A(u) texture(BuffD,(u)/RENDERSIZE.xy)

#define A(u) texture(BuffD,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffA,(u)/RENDERSIZE.xy)

#define A(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffA,(u)/RENDERSIZE.xy)
#define C(u) texture(image46,(u)/RENDERSIZE.xy)

#define A(u) texture(BuffB,(u)/RENDERSIZE.xy)
#define B(u) texture(BuffC,(u)/RENDERSIZE.xy)
*/
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;
    vec4  r = vec4(0);
    vec4  a = D(u);

    a.xy += _uvc*Succ*syn_Intensity;
    float z    = 4.;//kernel convolution size
    float blur = 4./z;
    blur *= 1.0 + _luminance( _edgeDetectSobel(syn_UserImage)*0.5*( 11.250+syn_BassLevel ));

    a += syn_Level*0.006;
    //a.r += (lum(_textureMedia(_uv))*4*PI)*(0.75*syn_BassLevel+0.35*syn_Intensity)*Media;
    //a.g += sin(lum(_loadUserImage())*4*PI)*(0.75*syn_BassLevel+0.35*syn_Intensity)*Media;

    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
        vec2  c = vec2(i,j)*blur; //c = c.yx*vec2(-1,1);
              c*= exp(-dot(c,c));
        vec4  a2= D(u+vec2(i,j));  
        vec4  b = a2-a;
        b *= 1.0 - ( _edgeDetectSobel(syn_UserImage)*0.35*( 1.0+syn_BassLevel ));
      

        r.x += length(b.xy-_uvc*Split);
        //r.xy += c*b.z;
        //r.z  += dot(c,b.xy           );//*a2.z;
        //r.w  += dot(c,b.yx*vec2(-1,1));//*a2.z;
    }}
    fragColor = r;
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;
    u += _uvc*PI*(Zoom+push_in*8. - push_out *8.);
    u += _uvc*PI*Stretch;
    u += Drift*PI;
    u += vec2(_noise(TIME)*0.5, 0.5*_noise(TIME))*syn_Presence;

    vec4  r = D(u);
    float z    = 4.;//kernel convolution size
    float blur = 3./z;
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec4 b = A(u+vec2(i,j));
      vec2 c = (-vec2(i,j))*blur;
           c*= exp(-dot(c,c))*.009;//*r.z;
      r.xy += +b.x*c              *1. //change to -1 to make it an atraction
              +b.w*c.yx*vec2(-1,1)*0.;
    }}

    fragColor = r;
	return fragColor; 
 } 


			//******** BuffC Code Begins ********

vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4 t = B(u);
    vec2 m = +t.xy
             +A(u).xy*(t.z-.5)*0.
             +t.z*vec2(0,.0)
             -T(u).x*t.xy*.0;
    float s = 0.;
    
    float z    = 4.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec2 c = (m+vec2(i,j))*1.;
      s += exp(-dot(c,c));
    }}
    if(s==0.){s = 1.;}
    s = 1./s;
    
    fragColor = vec4(m,s,0);
	return fragColor; 
 } 


			//******** BuffD Code Begins ********

vec4 renderPassD() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;
fragColor+= syn_HighLevel*0.125+syn_MidHighLevel*0.125;
    float tz = 0.;
    vec4 a = vec4(0);
    float z    = 6.;//kernel convolution size
    for(float i=-z; i<=z; ++i){
    for(float j=-z; j<=z; ++j){
      vec4 t = B(u+vec2(i,j)); t.z = 1.;
      vec4 m = C(u+vec2(i,j));
      //m *= 1.0 - ( _luminance(_edgeDetectSobel(syn_UserImage))*0.125*( 1.0+syn_BassLevel ));

      vec2 c = (m.xy-vec2(i,j))*1.;
      float z = t.z*exp(-dot(c,c));

      a.xy += z*m.xy;
      a.z  += z*m.z;

      tz   += z;
    }}
    if(tz==0.){tz = 1.;}
    a.xy /= tz;
    if(_mouse.z>0.)
    {
        vec2 m = 22.*(u-_mouse.xy)/RENDERSIZE.y;
        a += vec4(normalize(m),0,0)*exp(-dot(m,m))*.2;
    }
    if(FRAMECOUNT<=1 || Reset > 0.)
    {
        vec2 m = 4.*(u-RENDERSIZE.xy*.5)/RENDERSIZE.y;
        a = vec4(1,0,1,1)*exp(-dot(m,m));
    }
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    vec4 a = texture(BuffD,u);
    //a *= 1.0 +  (_edgeDetectSobel(syn_UserImage)*0.25*( 1.0+syn_BassLevel )),1, 1;

    //    a.r += sin(lum(_loadUserImage())*4*PI)*(0.75*syn_BassLevel+0.35*syn_Intensity)*Media;
    //a.g += sin(lum(_loadUserImage())*4*PI)*(0.75*syn_BassLevel+0.35*syn_Intensity)*Media;
    //a *= 1.0 -  _luminance( _edgeDetectSobelMedia()*0.25*( 1.0+syn_BassLevel )),1, 1;

    fragColor = a.z+a.z*sin(length(a.xy)+vec4(1,2,3,4)+0.);
    fragColor = +sin(a.x*2.+vec4(2,3,5,8)+0.)*(.25)
                +sin(a.y*2.+vec4(2,3,3,4)+0.)*.25*(1.)
                +.5*(1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}