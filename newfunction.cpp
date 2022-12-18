/*
ll pm(ll x,ll y)
{  ll res=1;   x%=mod1;
  while(y>0) 
  { if(y&1)res=(res*x)%mod1;
    y=y>>1;   x=(x*x)%mod1;}  return res;
}
*/
/* ----------------------------------------------------*/
/*
bool myc(pll p1,pll p2)
{
   if(p1.first!=p2.first)
    return p1.first<p2.first;
   else  return p1.second<p2.second;
}
vpll vp(ll a[],ll n)
{ vpll vt;
  lp(i,0,n)vt.eb(make_pair(a[i],i+1));
  sort(all(vt),myc); return vt;
}
*/
/*------------------------------------------------------*/
/*
vector<ll>primeser(ll n)
{ vector<ll>prime(n+1,1);
   prime[0]=0;prime[1]=0;
   lp(i,0,sqrt(n)+1)
   { if(prime[i])
     {for(ll j=i*i;j<=n;j+=i)
       prime[j]=0;}} return prime;
}
*/
/*------------------------------------------------------*/
/*
vll linearsieve (ll n) {
    vll prime;
    bool is_composite[n+2];
	fill (is_composite,is_composite+n,false);
	lp(i,2,n){
		if (!is_composite[i]) prime.eb(i);
		for (ll j=0;j<prime.size()&&i*prime[j]<n;++j) {
			is_composite[i * prime[j]] = true;
			if (i % prime[j] == 0) break;
		}
	}return prime;
}
*/
/*--------------------------------------------------------*/
/*
void subset(ll a[],ll n)
{
  lp(i,0,(1<<n)) {lp(j,0,n)
  { if(i&(1<<j))cout<<a[j]<<" ";}
    cout<<"\n";}
}
*/
/*--------------------------------------------------------*/
/*
bool is_prime(ll a)
{
  if(n==2||n==3)return true;
  if(n<=1||n%2==0||n%3==0)return false;
  lpa(i,5,sqrt(n)+1,6)if(n%i==0||n%(i + 2)==0)return false;
  return true;
}
*/
/*--------------------------------------------------------*/
/*
bool is_prime(ll a)
{
    if (a==2||a==3)return true;
    if (a<=1||a%2==0||a%3==0)return false;
    lpa(i,5,sqrt(a)+1,6)
    if (a%i==0||a%(i+2)== 0)return false;
    return true;
}
*/
/*--------------------------------------------------------*/
/*
ll minswap(ll a[],ll n)
{
   pll ap[n];
   lp(i,0,n){
    ap[i].first=a[i];ap[i].second=i;}
   sort(ap,ap+n);
   vector<bool>vis(n,false);
   ll ans=0;
   lp(i,0,n){
   if(vis[i]||ap[i].second==i)continue;
   ll cs=0;ll j=i;
    while(!vis[j]){
    vis[j]=1;j=ap[j].second;cs++;}
   if(cs)ans+=(cs-1);}
    return ans;
}
*/
/*----------------------------------------------------------*/
/*
ll merge(ll a[],ll temp[],ll le,ll mid,ll ri)
{ ll ic=0,i=le,j=mid,k=le;
  while((i<=mid-1)&&(j<=ri))
  { if(a[i]<=a[j])temp[k++]=a[i++];
    else{
    temp[k++]=a[j++];ic+=(mid-i);}}
  while(i<=mid-1)temp[k++]=a[i++];
  while(j<=ri)temp[k++]=a[j++];
  lp(i,le,ri+1)a[i]=temp[i];
   return ic;
}
ll mg(ll a[],ll temp[],ll le,ll ri)
{
   ll mid,ic=0;
   if(ri>le){
   mid=(ri+le)/2;
   ic=mg(a,temp,le,mid);
   ic+=mg(a,temp,mid+1,ri);
   ic+=merge(a,temp,le,mid+1,ri);}
   return ic;}  
ll countswap(ll a[],ll n)
{ ll temp[n];
  return mg(a,temp,0,n-1);
} 
*/
/*----------------------------------------------------------*/
/*
ll nopairsatiineq(ll l,ll r,vll &v,const ll diff){
         if(l==r)return 0;
         ll m=(l+r)/2;
         ll h=nopairsatiineq(l,m,v,diff)+nopairsatiineq(m+1,r,v,diff);
         vll vt;
         lp(j,m+1,r+1){
             ll k=v[j]+diff;
             ll p=upper_bound(v.begin()+l,v.begin()+m+1,k)-v.begin()-1;
             if(p<=m)h+=p-l+1;
         }
         int l1=l,r1=m,l2=m+1,r2=r;
         while(l1<=r1&&l2<=r2){
             if(v[l1]<=v[l2])vt.eb(v[l1++]);
             else vt.eb(v[l2++]); 
         }
         while(l1<=r1)vt.eb(v[l1++]);
         while(l2<=r2)vt.eb(v[l2++]);
         lp(i,l,r+1)v[i]=vt[i-l];
         return h;
     }    isme find karna hai number of i,j i.e. i<j and nums[i]-nums[j]<=diff
*/  
/*----------------------------------------------------------------*/
/*
ll kadane(ll a[],ll size)
{
    ll max_so_far =INM, max_ending_here = 0,start = 0, end = 0, s = 0;
    lp(i,0,size){
        max_ending_here += a[i];
        if (max_so_far < max_ending_here) {
            max_so_far = max_ending_here;
            start = s;
            end = i;
        }
        if (max_ending_here < 0) {
            max_ending_here = 0;
            s = i + 1;
        }
    }
    cout << "Maximum contiguous sum is " << max_so_far<< endl;
    cout << "Starting index " << start << endl<< "Ending index " << end << endl;
}
*/
/*-----------------------------------------------------------*/
/*
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        // http://xorshift.di.unimi.it/splitmix64.c
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};  iska use unordered map ke liye karenge sirf cp ke dauran hack se bachame ke liye
*/
/*---------------------------------------------------------------------*/
/*
void consttree(ll a[],ll seg[],ll l,ll h,ll p){
     if(l==h){seg[p]=a[l];return;}
     ll mid=(l+h)/2;
     consttree(a,seg,l,mid,2*p+1);
     consttree(a,seg,mid+1,h,2*p+2);
     seg[p]=min(seg[2*p+1],seg[2*p+2]);
 }
 ll rq(ll seg[],ll ql,ll qh,ll l,ll h,ll p){
       if(ql<=l&&qh>=h)return seg[p];
       if(ql>h||qh<l)return INF;
       ll mid=(l+h)/2;
       return min(rq(seg,ql,qh,l,mid,2*p+1),rq(seg,ql,qh,mid+1,h,2*p+2));
 }
 */
/*--------------------------------------------------------------------------*/
/*
ll ncr(ll n,ll r)
{
    if (r > n)return 0;
    ll inv[r + 1]={0},ans=1;
    inv[0]=1;
    if(r+1>=2)inv[1]=1;
    lp(i,2,r+1)inv[i]=mod2-(mod2/i)*inv[mod2%i]%mod2;
    lp(i,2,r+1)ans=((ans % mod2)*(inv[i]%mod2))%mod2;
    pl(i,n,n-r+1)ans=((ans%mod2)*(i % mod2))%mod2;
    return ans;
}
*/
/*--------------------------------------------------------------------------*/
/*
ll maxPoints(vector<vector<int>>& points) {
        ll n=points.size(),k=1;
        lp(i,0,n-1){
          map<pll,ll>m;
          lp(j,i+1,n){
              ll a=points[j][1]-points[i][1];
              ll b=points[j][0]-points[i][0];
              ll h=__gcd(a,b);
              a/=h;b/=h;
              m[{a,b}]++;
            }
            ll l=0;
            for(auto i:m)l=max(l,i.ss);
            k=max(k,l+1);
        }
        return k;
    }
*/
/*-------------------------------------------------------------------*/
/*
string nearestPalindromic(string n) {
        ll a=n.size(),res=0,mind=INF;
        vll v;
        v.eb((ll)(pow(10,a)+1));
        v.eb((ll)(pow(10,a-1)-1));
        ll k=((a+1)/2);
        ll h=stol(n.substr(0,k));
        vll vt={h,h+1,h-1};
        for(auto i:vt){
            string s=to_string(i);
            if(a&1)s.pop_back();
            reverse(all(s));
            string c=to_string(i)+s;
            v.eb(stol(c));
        }
        lp(i,0,5){
            if(v[i]!=stol(n)&&abs(v[i]-stol(n))<mind){
                mind=abs(v[i]-stol(n));
                res=v[i];
            }
            else if(abs(v[i]-stol(n))==mind)res=min(res,v[i]);
        }
        return to_string(res);
    }
*/
/*-----------------------------------------------------------------*/
/*
 ll countPairssuchthatproductdivisiblebyk(vll &nums,ll k) {
        ll n=nums.size(),ans=0;
        map<ll,ll>m;
        lp(i,0,n){
            ll h=__gcd(nums[i],k);
            for(auto j:m)if(j.ff*h%k==0)ans+=j.ss;
            m[h]++;
        }
        return ans;
    }
*/
/*------------------------------------------------------------------*/
/*void sleep(int delay)
 {   time_t t1,t2;
     t1=time(NULL);
     t2=t1+delay;
     while(t1<=t2)t1=time(NULL);
 }*/
 /*-------------------------------------------------------------------*/
 /*
 string DTB(ll n)
{
    //finding the binary form of the number and converting it to string. 
    string s = bitset<32> (n).to_string();
    return s;
}
*/
/*----------------------------------------------------------------------*/
/*
ll extendedEuclid(ll a,ll b,ll *x,ll *y){
    if(a>b)swap(a,b);
    if(a==0){
        *x=0;*y=1;
        return b;
    }
    ll x1,y1;
    ll g=extendedEuclid(b%a,a,&x1,&y1);
    *x=y1-(b/a)*x1;
    *y=x1;
    return g;
} 
here this x,y are used to find such that ax + by = gcd(a, b) 
Time Complexity: O(Log min(a, b))
*/
/*----------------------------------------------------------------------*/
/*
ll reduceB(ll a,string b)
{
    ll mod = 0;
    lp(i,0,b.size())mod=(mod*10+b[i]-'0')%a;
    return mod; // return modulo
}
 calculating mod of b with a to make b like 0 <= b < a
*/
/*----------------------------------------------------------------------*/
/*
ll steinsgcd(ll a,ll b)
{
    if (a == b)return a;
    if (a == 0)return b;
    if (b == 0)return a;
    if (~a & 1) 
    {  if (b & 1)return gcd(a >> 1, b);
        else return gcd(a >> 1, b >> 1) << 1;
    }
    if (~b & 1)return gcd(a, b >> 1);
    if(a<=b)swap(a,b);
    return gcd((a - b) >> 1, b);
}
O(N*N) where N is the number of bits in the larger number
*/
/*-----------------------------------------------------------------------*/
/*
ll PollardRho(ll n)
{
    srand (time(NULL));
    if (n==1) return n;
    if (n % 2 == 0) return 2;
    ll x = (rand()%(n-2))+2;
    ll y = x,d=1;
    ll c = (rand()%(n-1))+1;
    while (d==1)
    {
        x = (modular_pow(x, 2, n) + c + n)%n;
        y = (modular_pow(y, 2, n) + c + n)%n;
        y = (modular_pow(y, 2, n) + c + n)%n;
        d = __gcd(abs(x-y), n);
        if (d==n) return PollardRho(n);
    }
    return d;
} Given a positive integer n, and that it is composite, find a divisor of it. where prime
  can be very large
*/
/*----------------------------------------------------------------------*/
/*
void printDivisors(ll n)
{
    ll i;
    lp(i,1,sqrt(n)+1)if(n%i==0)cout<<i<<" ";
    if(i-(n/i)==1)i--;
    for(;i>=1;i--)if(n%i==0)cout<<n/i<<" ";
} Find all factors of a Natural Number in sorted order in auxilary space O(1)
*/
/*---------------------------------------------------------------------*/
/*
ll inv(ll a,ll m)
{
    ll m0=m,t,q,x0=0,x1=1;
    if (m==1)return 0;
    while (a > 1) {
        q = a / m;
        t = m;
        m=a%m,a=t;
        t=x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if(x1<0)x1+=m0;
    return x1;
}
ll findMinX(ll num[],ll rem[],ll k)
{
    ll prod = 1,res=0;
    lp(i,0,k)prod*=num[i];
    lp(i,0,k){
        ll pp=prod/num[i];
        res+=rem[i]*inv(pp,num[i])*pp;
    }
    return result % prod;
}  chinese remainder theorem
*/
/*----------------------------------------------------------------------------*/