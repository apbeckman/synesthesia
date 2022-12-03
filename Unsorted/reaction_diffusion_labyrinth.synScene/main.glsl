

			//******** BuffA Code Begins ********

#define A(u) texture(BuffA,(u)/RENDERSIZE.xy)
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 u = _xy;

    vec4 a = A(u+vec2(0,0));
    float b = 0.;
    
    //kernel convolution that blurs
    {
        float z = 4.;//kernel convolution size
        float t = z*2.+1.;
              t = 1./(t*t-1.);
        float blur = 1./z;
        for(float i = -z; i<=z;++i){
        for(float j = -z; j<=z;++j){
          float s = 0.;
          for(float i2 = -z; i2<=z;++i2){
          for(float j2 = -z; j2<=z;++j2){
          vec2 c = vec2(i2,j2)*blur;
          s += 1./exp(dot(c,c));
          }}
          if(s==0.)s = 1.;
          vec2 c = -vec2(i,j)*blur;
          float d = A(u+vec2(i,j)).x;
          b -= d/s/exp(dot(c,c))*9.;
        }}
    }
    //kernel convolution that blurs a little bit less
    {
        float z = 4.;//kernel convolution size
        float t = z*2.+1.;
              t = 1./(t*t-1.);
        float blur = 1.5/z;
        for(float i = -z; i<=z;++i){
        for(float j = -z; j<=z;++j){
          float s = 0.;
          for(float i2 = -z; i2<=z;++i2){
          for(float j2 = -z; j2<=z;++j2){
          vec2 c = vec2(i2,j2)*blur;
          s += 1./exp(dot(c,c));
          }}
          if(s==0.)s = 1.;
          vec2 c = -vec2(i,j)*blur;
          float d = A(u+vec2(i,j)).x;
          b += d/s/exp(dot(c,c))*9.8;
        }}
    }
    b = clamp(b,0.,1.);
    a = vec4(b,a.xyz);

    if(_mouse.z>0.)
    {
        vec2 m1 = 2.*(u-_mouse.xy)/RENDERSIZE.y;
        a *= 1.-+1./exp(pow(max(length(m1)-.1,0.),2.)*111.);
    }
    if(FRAMECOUNT==0)
    {
        vec2 m1 = (2.*u-RENDERSIZE.xy)/RENDERSIZE.y+vec2(.2,0);
        vec2 m2 = (2.*u-RENDERSIZE.xy)/RENDERSIZE.y-vec2(.2,0);
        a += +1./exp(pow(max(length(m1)-.2,0.),2.)*5555.);
        a += -1./exp(pow(max(length(m2)-.2,0.),2.)*5555.);
    }
    float keyA  = texture( iChannel1, vec2(65.5/256.,.25) ).x;
    if(keyA!=0.)a = texture( image5, u/RENDERSIZE.xy);
    fragColor = a;
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 u = fragCoord/RENDERSIZE.xy;
    fragColor = texture(BuffA,u);
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