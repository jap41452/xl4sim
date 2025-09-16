/* ===== DOM ===== */
const plateWidth   = document.getElementById('plateWidth');
const plateLength  = document.getElementById('plateLength');
const plateMaterial= document.getElementById('plateMaterial');
const coolantSel   = document.getElementById('coolant');
const flowLpm      = document.getElementById('flowLpm');
const inletT       = document.getElementById('inletT');

const gridCanvas   = document.getElementById('gridCanvas');
const sizeReadout  = document.getElementById('sizeReadout');
const tickReadout  = document.getElementById('tickReadout');

const maxCh    = document.getElementById('maxCh');
const minCh    = document.getElementById('minCh');
const chanB    = document.getElementById('chanB');
const minWall  = document.getElementById('minWall');
const chanH    = document.getElementById('chanH');
const plateThk = document.getElementById('plateThk');
const lamBC    = document.getElementById('lamBC');

const btnAddHS = document.getElementById('btnAddHS');
const btnClearHS = document.getElementById('btnClearHS');
const hsList   = document.getElementById('hsList');
const totalHeat= document.getElementById('totalHeat');

const btnSolve = document.getElementById('btnSolve');
const btnContours = document.getElementById('btnContours');
const chart    = document.getElementById('chart');
const chartLegend = document.getElementById('chartLegend');
const solveNote = document.getElementById('solveNote');

const btnCSV   = document.getElementById('btnCSV');
const btnJSON  = document.getElementById('btnJSON');

const endView  = document.getElementById('endView');
const endSummary = document.getElementById('endSummary');
const endNote    = document.getElementById('endNote');

const statusBody = document.getElementById('statusBody');

const contourCanvas = document.getElementById('contourCanvas');
const contourSummary = document.getElementById('contourSummary');
const contourNote = document.getElementById('contourNote');

/* ===== State ===== */
let heatSources = []; // {id,x,y,w,h,W}
let drawing=false, dragStart=null, pendingAdd=false;
let lastSweep=null; // {Ns, TmaxC, dPkPa, Re, Nu, h, bcUsed, Tfields:{[N]:{nx,ny,L,W,T,Tmin,Tmax,bc}}}
let selectionMode=false;
let chartPts=[]; // {x,y,i}

/* ===== Helpers ===== */
function parseNum(el,fallback){ const v=parseFloat(el.value); return (isFinite(v)?v:fallback); }
function niceStep(maxVal,targetTicks=9){
  const raw=maxVal/targetTicks, p10=Math.pow(10,Math.floor(Math.log10(raw)));
  const n=raw/p10; let m; if(n<1.5)m=1; else if(n<3.5)m=2; else if(n<7.5)m=5; else m=10; return m*p10;
}
const mm2m=function(v){ return v/1000; };
const clamp=function(v,a,b){ return Math.max(a,Math.min(b,v)); };

/* ===== Materials & Fluids ===== */
function plateK(material){ switch(material){ case 'Copper': return 390; case 'Stainless Steel': return 16; case 'FRP': return 0.3; default: return 180; } }
function fluidProps(name){
  if(name==='Glycol 50/50 mix'){ const rho=1040, mu=0.0049, k=0.37, cp=3300; return {rho:rho,mu:mu,k:k,cp:cp,Pr:(mu*cp/k)}; }
  const rho=997, mu=0.00089, k=0.6, cp=4180; return {rho:rho,mu:mu,k:k,cp:cp,Pr:(mu*cp/k)};
}

/* ===== Grid & transforms (Plate Grid) ===== */
let x0=0,y0=0,plotW=0,plotH=0,Lmm=200,Wmm=200;
function xToPx(x){ return x0 + (x / Lmm) * plotW; }
function yToPx(y){ return (y0 + plotH) - (y / Wmm) * plotH; }
function pxToX(px){ return ((px - x0) / plotW) * Lmm; }
function pxToY(py){ return (((y0 + plotH) - py) / plotH) * Wmm; }

function drawGrid(){
  Wmm=parseNum(plateWidth,200); Lmm=parseNum(plateLength,200);
  const panel=gridCanvas.parentElement.parentElement;
  const pad=parseFloat(getComputedStyle(panel).paddingLeft)+parseFloat(getComputedStyle(panel).paddingRight);
  const innerW=panel.clientWidth - pad;

  const m={left:36,right:8,top:10,bottom:22};
  plotW=Math.max(10, innerW - (m.left+m.right));
  plotH=Math.max(10, plotW * (Wmm/Lmm));

  gridCanvas.width=Math.round(m.left+plotW+m.right);
  gridCanvas.height=Math.round(m.top+plotH+m.bottom);

  const ctx=gridCanvas.getContext('2d');
  ctx.clearRect(0,0,gridCanvas.width,gridCanvas.height);
  ctx.fillStyle='#fff';ctx.fillRect(0,0,gridCanvas.width,gridCanvas.height);
  ctx.strokeStyle='#d9d9d9';ctx.strokeRect(0.5,0.5,gridCanvas.width-1,gridCanvas.height-1);

  x0=m.left;y0=m.top; const x1=x0+plotW, y1=y0+plotH;
  const majorX=niceStep(Lmm,9), majorY=niceStep(Wmm,9), minorDiv=5;

  ctx.save();ctx.lineWidth=1;ctx.strokeStyle='#eeeeee';ctx.beginPath();
  for(let x=0;x<=Lmm+1e-6;x+=majorX/minorDiv){const px=Math.round(xToPx(x))+0.5;ctx.moveTo(px,y0);ctx.lineTo(px,y1);}
  for(let y=0;y<=Wmm+1e-6;y+=majorY/minorDiv){const py=Math.round(yToPx(y))+0.5;ctx.moveTo(x0,py);ctx.lineTo(x1,py);} ctx.stroke();ctx.restore();

  ctx.save();ctx.lineWidth=1;ctx.strokeStyle='#cfd4dc';ctx.beginPath();
  for(let x=0;x<=Lmm+1e-6;x+=majorX){const px=Math.round(xToPx(x))+0.5;ctx.moveTo(px,y0);ctx.lineTo(px,y1);}
  for(let y=0;y<=Wmm+1e-6;y+=majorY){const py=Math.round(yToPx(y))+0.5;ctx.moveTo(x0,py);ctx.lineTo(x1,py);} ctx.stroke();ctx.restore();

  ctx.save();ctx.lineWidth=2;ctx.strokeStyle='#0b73c8';
  ctx.strokeRect(Math.round(x0)+0.5, Math.round(y0)+0.5, Math.round(plotW), Math.round(plotH));
  ctx.restore();

  ctx.save();ctx.fillStyle='#111';ctx.font='11px system-ui,-apple-system,Segoe UI,Roboto,sans-serif';
  ctx.textAlign='center';ctx.textBaseline='top';
  for(let x=0;x<=Lmm+1e-6;x+=majorX){ctx.fillText(x.toFixed(0), xToPx(x), y1+4);}
  ctx.textAlign='right';ctx.textBaseline='middle';
  for(let y=0;y<=Wmm+1e-6;y+=majorY){ctx.fillText(y.toFixed(0), x0-6, yToPx(y));}
  ctx.fillStyle='#666';ctx.textAlign='center';ctx.fillText('x (mm)', (x0+x1)/2, y1+16);
  ctx.save();ctx.translate(x0-26, (y0+y1)/2);ctx.rotate(-Math.PI/2);ctx.fillText('y (mm)', 0, 0);ctx.restore();ctx.restore();

  drawHeatOverlays(ctx);
  sizeReadout.textContent='Plate: ' + Wmm + ' mm (W) × ' + Lmm + ' mm (L)';
  tickReadout.textContent='Major: Δx=' + majorX + ' mm, Δy=' + majorY + ' mm';
  drawEndView();
}

/* ===== Heat sources UI/overlay ===== */
const HS_COLORS=['#e11d48','#0ea5e9','#22c55e','#f59e0b','#8b5cf6'];
function drawHeatOverlays(ctx){
  heatSources.forEach(function(hs,idx){
    const color=HS_COLORS[idx%HS_COLORS.length];
    const x=xToPx(hs.x), y=yToPx(hs.y+hs.h);
    const w=xToPx(hs.x+hs.w)-xToPx(hs.x);
    const h=yToPx(hs.y)-yToPx(hs.y+hs.h);
    ctx.save();
    ctx.fillStyle=hexToRgba(color,0.18); ctx.strokeStyle=color; ctx.lineWidth=2;
    ctx.fillRect(x,y,w,h); ctx.strokeRect(x+0.5,y+0.5,w-1,h-1);
    ctx.fillStyle=color; ctx.fillRect(x+4,y+4,24,16);
    ctx.fillStyle='#fff'; ctx.font='11px system-ui,sans-serif'; ctx.textBaseline='top'; ctx.textAlign='left'; ctx.fillText(String(idx+1), x+9, y+5);
    ctx.fillStyle='#111'; ctx.textBaseline='bottom'; ctx.textAlign='left'; ctx.fillText(hs.W.toFixed(0) + ' W', x+6, y+h-6);
    ctx.restore();
  });
}
function hexToRgba(hex,a){const m=/^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex); if(!m) return 'rgba(0,0,0,'+a+')';
  const r=parseInt(m[1],16), g=parseInt(m[2],16), b=parseInt(m[3],16); return 'rgba('+r+','+g+','+b+','+a+')'; }

btnAddHS.addEventListener('click', function(){ if(heatSources.length>=5){alert('Limit reached: up to 5 heat sources.'); return;}
  pendingAdd=true; btnAddHS.textContent='Adding… (drag on grid)'; btnAddHS.disabled=true; });
btnClearHS.addEventListener('click', function(){ heatSources=[]; refreshHSList(); updateTotalHeat(); drawGrid(); });

gridCanvas.addEventListener('mousedown', function(e){ if(!pendingAdd) return; const r=gridCanvas.getBoundingClientRect();
  dragStart={x:e.clientX-r.left,y:e.clientY-r.top}; drawing=true; });
gridCanvas.addEventListener('mousemove', function(e){ if(!drawing) return; drawGrid();
  const r=gridCanvas.getBoundingClientRect(); const x2=e.clientX-r.left,y2=e.clientY-r.top; const x1=dragStart.x,y1=dragStart.y;
  const ctx=gridCanvas.getContext('2d'); ctx.save(); ctx.strokeStyle='#0b73c8'; ctx.setLineDash([6,4]); ctx.lineWidth=1.5;
  ctx.strokeRect(Math.min(x1,x2),Math.min(y1,y2),Math.abs(x2-x1),Math.abs(y2-y1)); ctx.restore(); });
gridCanvas.addEventListener('mouseup', function(e){
  if(!drawing) return; drawing=false; const r=gridCanvas.getBoundingClientRect(); const x2=e.clientX-r.left,y2=e.clientY-r.top;
  const x1=dragStart.x,y1=dragStart.y; dragStart=null;
  const minPxX=Math.max(Math.min(x1,x2),x0), maxPxX=Math.min(Math.max(x1,x2),x0+plotW);
  const minPxY=Math.max(Math.min(y1,y2),y0), maxPxY=Math.min(Math.max(y1,y2),y0+plotH);
  if(maxPxX-minPxX<3 || maxPxY-minPxY<3){ drawGrid(); finishAddMode(); return; }
  const x_mm=Math.max(0, pxToX(minPxX));
  const y_mm_bottom=Math.max(0, pxToY(maxPxY));
  const w_mm=Math.min(Lmm, pxToX(maxPxX)) - x_mm;
  const h_mm=Math.min(Wmm, pxToY(minPxY)) - y_mm_bottom;
  let W=prompt('Enter heat source power (W):','250'); if(W===null){ drawGrid(); finishAddMode(); return; }
  W=parseFloat(W); if(!isFinite(W)||W<=0){ alert('Please enter a positive number.'); drawGrid(); finishAddMode(); return; }
  heatSources.push({id:rid(), x:x_mm, y:y_mm_bottom, w:w_mm, h:h_mm, W:W});
  refreshHSList(); updateTotalHeat(); drawGrid(); finishAddMode();
});
function finishAddMode(){ pendingAdd=false; btnAddHS.textContent='Add'; btnAddHS.disabled=false; }
function rid(){ return 'hs_'+Math.random().toString(36).slice(2,9); }
function refreshHSList(){ hsList.innerHTML='';
  heatSources.forEach(function(hs,idx){
    const row=document.createElement('div'); row.className='hs-item';
    const name=document.createElement('div'); name.textContent='#'+(idx+1)+' – '+hs.W.toFixed(0)+' W';
    const dims=document.createElement('span'); dims.className='pill'; dims.textContent=hs.w.toFixed(1)+'×'+hs.h.toFixed(1)+' mm';
    const del=document.createElement('button'); del.className='del'; del.textContent='×'; del.title='Remove';
    del.onclick=function(){ heatSources.splice(idx,1); refreshHSList(); updateTotalHeat(); drawGrid(); };
    row.appendChild(name); row.appendChild(dims); row.appendChild(del); hsList.appendChild(row);
  });
}
function updateTotalHeat(){ const sum=heatSources.reduce(function(s,h){return s+h.W;},0); totalHeat.value=sum.toFixed(1); }

/* ===== Channel counts & end view (fixed walls) ===== */
function updateChannelCounts(){ const W=parseNum(plateWidth,200); maxCh.value=Math.floor(W/3); minCh.value=Math.floor(W/8); }
function channelWidthMM(Wmm,N,tmm){ return (Wmm - (N+1)*tmm)/N; }

function drawEndView(){
  const Wmm_ = parseNum(plateWidth,200);
  const tmm  = Math.max(0, parseFloat(minWall.value)||0);
  const bmm  = Math.max(0, parseFloat(chanB.value)||0);
  const tplate = Math.max(0.1, parseFloat(plateThk.value)||2.0);

  let Nmin = Math.max(1, Math.floor(Wmm_/8));
  let N = Nmin, wmm = channelWidthMM(Wmm_, N, tmm);
  while (N>1 && wmm<=0){ N--; wmm=channelWidthMM(Wmm_, N, tmm); }
  const feasible = wmm>0;
  chanH.value = feasible ? wmm.toFixed(2) : '—';

  const panel=endView.parentElement; const pad=parseFloat(getComputedStyle(panel).paddingLeft)+parseFloat(getComputedStyle(panel).paddingRight);
  const innerW=panel.clientWidth - pad;
  const m={left:36,right:8,top:12,bottom:32}; const plateHpx=70;
  endView.width=Math.round(m.left+(innerW-(m.left+m.right))+m.right);
  endView.height=Math.round(m.top+plateHpx+18+m.bottom);
  const ctx=endView.getContext('2d');
  ctx.clearRect(0,0,endView.width,endView.height);
  ctx.fillStyle='#fff'; ctx.fillRect(0,0,endView.width,endView.height);
  ctx.strokeStyle='#d9d9d9'; ctx.strokeRect(0.5,0.5,endView.width-1,endView.height-1);
  const plotWpx=endView.width-(m.left+m.right);
  const X=function(mm){ return m.left + (mm/Wmm_)*plotWpx; };
  const yTop=m.top, yBot=yTop+plateHpx;

  ctx.fillStyle='#f8fafc'; ctx.fillRect(X(0), yTop, X(Wmm_)-X(0), plateHpx);
  ctx.strokeStyle='#0b73c8'; ctx.lineWidth=2; ctx.strokeRect(Math.round(X(0))+0.5, Math.round(yTop)+0.5, Math.round(X(Wmm_)-X(0)), Math.round(plateHpx)-1);

  if (feasible){
    const depthRatio=Math.min(1, bmm/tplate); const bpx=plateHpx*depthRatio;
    let x=tmm; for(let k=0;k<N;k++){ const x0=x, x1=x+wmm;
      ctx.fillStyle='rgba(56,189,248,0.25)'; ctx.strokeStyle='#0ea5e9'; ctx.lineWidth=1.5;
      ctx.fillRect(X(x0), yTop, X(x1)-X(x0), bpx);
      ctx.strokeRect(X(x0)+0.5, yTop+0.5, X(x1)-X(x0)-1, bpx-1);
      x = x1 + tmm;
    }
  }

  const major=niceStep(Wmm_,9);
  ctx.save(); ctx.strokeStyle='#cfd4dc'; ctx.lineWidth=1; ctx.beginPath();
  for(let x=0;x<=Wmm_+1e-6;x+=major){ const px=Math.round(X(x))+0.5; ctx.moveTo(px,yBot+4); ctx.lineTo(px,yBot+10); }
  ctx.stroke();
  ctx.fillStyle='#111'; ctx.font='11px system-ui,sans-serif'; ctx.textAlign='center'; ctx.textBaseline='top';
  for(let x=0;x<=Wmm_+1e-6;x+=major) ctx.fillText(String(x.toFixed(0)), X(x), yBot+12);
  ctx.fillStyle='#666'; ctx.fillText('width (mm)', X(Wmm_/2), yBot+26); ctx.restore();

  endSummary.textContent = feasible
    ? 'N(min)='+Nmin+', used='+N+' • t='+tmm.toFixed(2)+' mm • w='+wmm.toFixed(2)+' mm • depth b='+bmm.toFixed(2)+' mm • plate='+tplate.toFixed(2)+' mm'
    : 'Not feasible: (W − (N+1)t)/N ≤ 0 for N='+Nmin;
  endNote.textContent = (feasible && N < Nmin) ? 'Reduced N to fit fixed walls' : '';
}

/* ===== Nu correlations & hydraulics ===== */
function Nu_lam_rect_H2_allIsoflux(w,b){
  const alpha = Math.min(w,b)/Math.max(w,b);
  const beta  = (1 - alpha) / (1 + alpha);
  const br = 1 - 2.0421*beta + 3.0853*beta*beta - 2.4765*beta*beta*beta + 1.0578*Math.pow(beta,4) - 0.1861*Math.pow(beta,5);
  return 8.235 * br;
}
function Nu_lam_rect_T_allIso(w,b){ return 0.915 * Nu_lam_rect_H2_allIsoflux(w,b); }
function Nu_lam_rect_H2_topOnly(w,b){ return 1.18 * Nu_lam_rect_H2_allIsoflux(w,b); }

function frictionFactor(Re){ if(Re<2300) return 64/Math.max(Re,1e-12); const inv=-1.8*Math.log10(6.9/Math.max(Re,1e-12)); return 1/(inv*inv); }
function Nu_turb_DB(Re,Pr){ return 0.023*Math.pow(Re,0.8)*Math.pow(Pr,0.4); }

/* ===== Plate source & masks ===== */
function buildQppField(nx,ny,L,W){
  const qpp=Array.from({length:ny},function(){return Array(nx).fill(0);});
  if(heatSources.length===0) return qpp;
  const dx=L/nx, dy=W/ny;
  for(let hsi=0; hsi<heatSources.length; hsi++){
    const hs = heatSources[hsi];
    const Ax=mm2m(hs.w)*mm2m(hs.h); if(Ax<=0) continue;
    const qhs=hs.W/Ax;
    for(let j=0;j<ny;j++){ const yc=(j+0.5)*dy; if(yc<mm2m(hs.y)||yc>mm2m(hs.y+hs.h)) continue;
      for(let i=0;i<nx;i++){ const xc=(i+0.5)*dx; if(xc<mm2m(hs.x)||xc>mm2m(hs.x+hs.w)) continue; qpp[j][i]+=qhs; }
    }
  }
  return qpp;
}
function buildMasksFixedWall(N, nx, ny, W, w, t){
  const masks=Array.from({length:N},function(){return Array.from({length:ny},function(){return Array(nx).fill(false);});});
  const dy=W/ny; const starts=[]; let x=t; for(let k=0;k<N;k++){ starts.push([x,x+w]); x+=w+t; }
  for(let k=0;k<N;k++){ const pair=starts[k]; const x0=pair[0], x1=pair[1];
    for(let j=0;j<ny;j++){ const yc=(j+0.5)*dy; if(yc>=x0 && yc<=x1){ for(let i=0;i<nx;i++) masks[k][j][i]=true; } }
  }
  return masks;
}

/* ===== Linear solver (matrix-free CG) ===== */
const dot=function(a,b){let s=0;for(let i=0;i<a.length;i++) s+=a[i]*b[i];return s;};
const idx=function(i,j,nx){ return j*nx + i; };
function cgSolve(applyA,b,x0,maxIt,tol){
  if(typeof maxIt==='undefined') maxIt=4000;
  if(typeof tol==='undefined') tol=1e-6;
  const n=b.length; const x=x0.slice(); const r=new Array(n),p=new Array(n),Ap=new Array(n);
  applyA(x,Ap); for(let i=0;i<n;i++){ r[i]=b[i]-Ap[i]; p[i]=r[i]; }
  let rs=dot(r,r), bnorm=Math.max(1e-12, Math.sqrt(dot(b,b)));
  for(let it=0; it<maxIt; it++){
    applyA(p,Ap); const alpha = rs/Math.max(1e-30, dot(p,Ap));
    for(let i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }
    const rsn=dot(r,r); if(Math.sqrt(rsn)/bnorm < tol) break;
    const beta=rsn/Math.max(1e-30,rs); for(let i=0;i<n;i++) p[i]=r[i]+beta*p[i]; rs=rsn;
  }
  return x;
}
function buildOperator(nx,ny,L,W,ks,t,h,masks,Tb){
  const dx=L/nx, dy=W/ny, Acell=dx*dy, kt=ks*t;
  const aE=kt*(dy/dx), aW=aE, aN=kt*(dx/dy), aS=aN;

  const hA=Array.from({length:ny},function(){return Array(nx).fill(0);});
  const hTbA=Array.from({length:ny},function(){return Array(nx).fill(0);});
  for(let k=0;k<masks.length;k++){
    for(let i=0;i<nx;i++){
      const Tb_i=Tb[k][i];
      for(let j=0;j<ny;j++){
        if(masks[k][j][i]){ hA[j][i]+= h*Acell; hTbA[j][i]+= h*Acell*Tb_i; }
      }
    }
  }
  function assembleB(qpp){
    const b=new Array(nx*ny);
    for(let j=0;j<ny;j++) for(let i=0;i<nx;i++) b[idx(i,j,nx)] = qpp[j][i]*Acell + hTbA[j][i];
    return b;
  }
  function applyA(x,y){
    for(let j=0;j<ny;j++){
      for(let i=0;i<nx;i++){
        const p=idx(i,j,nx);
        let diag=hA[j][i], acc=0;
        if(i<nx-1){ acc+=aE * x[idx(i+1,j,nx)]; diag+=aE; }
        if(i>0   ){ acc+=aW * x[idx(i-1,j,nx)]; diag+=aW; }
        if(j>0   ){ acc+=aN * x[idx(i,  j-1,nx)]; diag+=aN; }
        if(j<ny-1){ acc+=aS * x[idx(i,  j+1,nx)]; diag+=aS; }
        y[p] = diag * x[p] - acc;
      }
    }
  }
  return { assembleB: assembleB, applyA: applyA };
}
function updateTbFromPlate(nx,ny,L,W,h,masks,T,TbOld,mDotPer,cp,relax){
  if(typeof relax==='undefined') relax=0.5;
  const dx=L/nx, dy=W/ny, Acell=dx*dy, K=masks.length;
  const TbNew=Array.from({length:K},function(){return Array(nx+1).fill(TbOld[0][0]);});
  for(let k=0;k<K;k++){
    TbNew[k][0]=TbOld[k][0];
    for(let i=0;i<nx;i++){
      let qx=0; for(let j=0;j<ny;j++) if(masks[k][j][i]) qx += h*(T[j][i]-TbOld[k][i])*Acell;
      const dT=(mDotPer>0 ? qx/(mDotPer*cp) : 0);
      const cand = TbNew[k][i] + dT;
      TbNew[k][i+1] = relax*cand + (1-relax)*TbOld[k][i+1];
    }
  }
  return TbNew;
}
function initialGuess(nx,ny,T0){ const a=new Array(nx*ny); for(let j=0;j<ny;j++) for(let i=0;i<nx;i++) a[idx(i,j,nx)]=T0; return a; }
function reshape(flat,nx,ny){const A=Array.from({length:ny},function(){return Array(nx);}); for(let j=0;j<ny;j++) for(let i=0;i<nx;i++) A[j][i]=flat[idx(i,j,nx)]; return A;}
function max2D(A){ let m=-1e9; for(let r=0;r<A.length;r++){ for(let c=0;c<A[r].length;c++){ if(A[r][c]>m) m=A[r][c]; } } return m; }
function min2D(A){ let m=1e9; for(let r=0;r<A.length;r++){ for(let c=0;c<A[r].length;c++){ if(A[r][c]<m) m=A[r][c]; } } return m; }

/* ===== Solve & sweep ===== */
btnSolve.addEventListener('click', runSolve);
btnContours.addEventListener('click', enterContourPickMode);

function runSolve(){
  btnSolve.disabled=true; btnContours.disabled=true; solveNote.textContent='Solving…';
  setTimeout(function(){ lastSweep=sweepFiveCases(); solveNote.textContent='Solved N = ['+lastSweep.Ns.join(', ')+']'; btnSolve.disabled=false; btnContours.disabled=false; drawChart(lastSweep); renderStatus(lastSweep); if(selectionMode) hintPick(); }, 10);
}

function sweepFiveCases(){
  const Wmm=parseNum(plateWidth,200), Lmm=parseNum(plateLength,200);
  const W=mm2m(Wmm), L=mm2m(Lmm);
  const tmm=Math.max(0, parseFloat(minWall.value)||0);
  const bmm=Math.max(0, parseFloat(chanB.value)||0);
  const b=mm2m(bmm);
  const t_plate=mm2m(Math.max(0.2, parseFloat(plateThk.value)||2.0));
  const ks=plateK(plateMaterial.value);
  const fluid=fluidProps(coolantSel.value);
  const Q_tot=Math.max(0, parseFloat(flowLpm.value)||0)/1000/60;
  const Tin=parseFloat(inletT.value)||25;
  const nx=120, ny=Math.max(12, Math.round(nx*(W/L)));
  const qpp=buildQppField(nx,ny,L,W);

  let Nmin=Math.max(1, Math.floor(Wmm/8));
  let Nmax=Math.max(1, Math.floor(Wmm/3)); if(Nmax<Nmin) { const tmp=Nmax; Nmax=Nmin; Nmin=tmp; }
  const span=Nmax-Nmin, picks=[0,1,2,3,4].map(function(k){ return Math.round(Nmin+(span*k/4)); });
  let Ns=Array.from(new Set(picks)).filter(function(n){return n>=1;});
  Ns = Ns.map(function(n){ let wmm2 = channelWidthMM(Wmm,n,tmm); while(n>1 && wmm2<=0){ n--; wmm2=channelWidthMM(Wmm,n,tmm); } return Math.max(1,n); });
  Ns = Array.from(new Set(Ns)).sort(function(a,b){return a-b;});

  const TmaxC=[], dPkPa=[], ReA=[], NuA=[], hA=[], bcUsed=[], Tfields={};

  for(let idxN=0; idxN<Ns.length; idxN++){
    const N = Ns[idxN];
    const wmm = channelWidthMM(Wmm, N, tmm);
    const w   = mm2m(wmm);
    const A   = w*b;
    const Dh  = (A>0 ? (2*w*b)/(w+b) : 1e-6);
    const Qch = (N>0? Q_tot/N : 0);
    const v   = (A>0? Qch/A : 0);

    const Re  = fluid.rho * v * Dh / Math.max(fluid.mu,1e-12);
    let Nu, bcTag;
    if (Re < 2300){
      const bc = lamBC.value;
      if (bc==='T'){ Nu = Nu_lam_rect_T_allIso(w,b); bcTag='T'; }
      else if (bc==='H2-1'){ Nu = Nu_lam_rect_H2_topOnly(w,b); bcTag='H2-1'; }
      else { Nu = Nu_lam_rect_H2_allIsoflux(w,b); bcTag='H2-4'; }
    } else { Nu = Nu_turb_DB(Re, fluid.Pr); bcTag='Turb'; }
    const h   = Nu * fluid.k / Math.max(Dh,1e-9);

    const f   = frictionFactor(Re);
    const dP  = (f*(L/Dh))*(fluid.rho*v*v/2); // Pa
    dPkPa.push(dP/1000); ReA.push(Re); NuA.push(Nu); hA.push(h); bcUsed.push(bcTag);

    const masks = buildMasksFixedWall(N, nx, ny, W, w, mm2m(tmm));
    let Tb = Array.from({length:N},function(){return Array(nx+1).fill(Tin);});
    let T  = Array.from({length:ny},function(){return Array(nx).fill(Tin);});

    var op0 = buildOperator(nx,ny,L,W,ks,t_plate,h,masks,Tb);
    var assembleB = op0.assembleB;
    var applyA = op0.applyA;
    let bvec = assembleB(qpp), x = cgSolve(applyA,bvec, initialGuess(nx,ny,Tin), 5000, 1e-6);
    T = reshape(x,nx,ny);
    for(let it=0; it<2; it++){
      Tb = updateTbFromPlate(nx,ny,L,W,h,masks,T,Tb, fluid.rho*Qch, fluid.cp, 0.5);
      var opi = buildOperator(nx,ny,L,W,ks,t_plate,h,masks,Tb);
      assembleB = opi.assembleB;
      applyA   = opi.applyA;
      bvec = assembleB(qpp); x = cgSolve(applyA,bvec,x, 4000, 1e-6);
      T = reshape(x,nx,ny);
    }

    const Tmax = max2D(T), Tmin = min2D(T);
    TmaxC.push(Tmax);
    Tfields[N] = { nx: nx, ny: ny, L: L, W: W, T: T, Tmin: Tmin, Tmax: Tmax, N: N, bc: bcTag };
  }

  const w_ui = channelWidthMM(Wmm, Ns[0], tmm);
  chanH.value = isFinite(w_ui) && w_ui>0 ? w_ui.toFixed(2) : '—';

  return { Ns: Ns, TmaxC: TmaxC, dPkPa: dPkPa, Re: ReA, Nu: NuA, h: hA, bcUsed: bcUsed, Tfields: Tfields };
}

/* ===== Chart + pick mode ===== */
function drawChart(sweep){
  if(!sweep) return; chart.style.display='block'; chartLegend.style.display='flex';
  const ctx=chart.getContext('2d'); const Wc=chart.width, Hc=chart.height; ctx.clearRect(0,0,Wc,Hc);
  const pad={l:34,r:34,t:10,b:24}, pl=pad.l, pr=Wc-pad.r, pt=pad.t, pb=Hc-pad.b;
  const Ns=sweep.Ns, Tmax=sweep.TmaxC, DP=sweep.dPkPa;
  const xMin=Math.min.apply(null,Ns), xMax=Math.max.apply(null,Ns);
  const yLmin=Math.min.apply(null,Tmax), yLmax=Math.max.apply(null,Tmax);
  const yRmin=Math.min.apply(null,DP),   yRmax=Math.max.apply(null,DP);
  const x=function(N){return pl+(pr-pl)*((N-xMin)/Math.max(1e-9,(xMax-xMin)));};
  const yL=function(v){return pb-(pb-pt)*((v-yLmin)/Math.max(1e-9,(yLmax-yLmin)));};
  const yR=function(v){return pb-(pb-pt)*((v-yRmin)/Math.max(1e-9,(yRmax-yRmin)));};

  // grid
  ctx.strokeStyle='#e5e7eb'; ctx.lineWidth=1; ctx.beginPath();
  Ns.forEach(function(N){ ctx.moveTo(x(N)+0.5,pt); ctx.lineTo(x(N)+0.5,pb); }); ctx.stroke();
  ctx.strokeStyle='#111'; ctx.lineWidth=1.5; ctx.strokeRect(pl+0.5,pt+0.5,(pr-pl)-1,(pb-pt)-1);

  // left series: Tmax
  ctx.strokeStyle='#0b73c8'; ctx.lineWidth=2; ctx.beginPath();
  Ns.forEach(function(N,i){ const X=x(N), Y=yL(Tmax[i]); if(i===0) ctx.moveTo(X,Y); else ctx.lineTo(X,Y); }); ctx.stroke();
  ctx.fillStyle='#0b73c8'; Ns.forEach(function(N,i){ ctx.beginPath(); ctx.arc(x(N),yL(Tmax[i]),3,0,2*Math.PI); ctx.fill(); });

  // right series: DP
  ctx.strokeStyle='#e11d48'; ctx.lineWidth=2; ctx.beginPath();
  Ns.forEach(function(N,i){ const X=x(N), Y=yR(DP[i]); if(i===0) ctx.moveTo(X,Y); else ctx.lineTo(X,Y); }); ctx.stroke();
  ctx.fillStyle='#e11d48'; Ns.forEach(function(N,i){ ctx.beginPath(); ctx.rect(x(N)-3,yR(DP[i])-3,6,6); ctx.fill(); });

  // axes labels
  ctx.fillStyle='#111'; ctx.font='12px system-ui,sans-serif'; ctx.textAlign='center'; ctx.textBaseline='top';
  Ns.forEach(function(N){ctx.fillText(String(N), x(N), pb+4);}); ctx.fillStyle='#555'; ctx.fillText('N (channels)', (pl+pr)/2, Hc-14);
  ctx.fillStyle='#0b73c8'; ctx.textAlign='right'; ctx.textBaseline='middle';
  for(let i=0;i<=4;i++){ const v=yLmin+(i/4)*(yLmax-yLmin); ctx.fillText(v.toFixed(1), pl-6, yL(v)); }
  ctx.fillStyle='#e11d48'; ctx.textAlign='left'; ctx.textBaseline='middle';
  for(let i=0;i<=4;i++){ const v=yRmin+(i/4)*(yRmax-yRmin); ctx.fillText(v.toFixed(2), pr+6, yR(v)); }

  // cache left-series points for picking
  chartPts = Ns.map(function(N, i){
    return { x: x(N), y: yL(Tmax[i]), i: i };
  });
}

chart.addEventListener('click', function(e){
  if(!selectionMode || !lastSweep) return;
  const rect=chart.getBoundingClientRect();
  const mx=e.clientX-rect.left, my=e.clientY-rect.top;
  let best=-1, bd=1e9;
  chartPts.forEach(function(p){ const d=Math.hypot(mx-p.x,my-p.y); if(d<bd){ bd=d; best=p.i; } });
  if(best>=0 && bd<=10){
    const N=lastSweep.Ns[best];
    const fld=lastSweep.Tfields[N];
    if(fld) { renderContourPanel(fld); selectionMode=false; chart.classList.remove('pick'); contourNote.textContent=''; }
  } else {
    contourNote.textContent='Click directly on a Max-T point (blue dots).';
  }
});

function enterContourPickMode(){
  if(!lastSweep){ runSolve(); }
  selectionMode=true;
  hintPick();
}
function hintPick(){
  chart.classList.add('pick');
  contourNote.textContent='Pick a Max-T point on the chart to render contours.';
}

/* ===== Status table ===== */
function renderStatus(sweep){
  statusBody.innerHTML='';
  if(!sweep) return;
  for(let i=0;i<sweep.Ns.length;i++){
    const tr=document.createElement('tr');
    const cells=[ sweep.Ns[i], sweep.bcUsed[i], sweep.Re[i].toFixed(0), sweep.Nu[i].toFixed(2), sweep.h[i].toFixed(0), sweep.dPkPa[i].toFixed(2), sweep.TmaxC[i].toFixed(2) ];
    cells.forEach(function(txt,ci){ const td=document.createElement('td'); td.textContent=txt; if(ci===0) td.style.textAlign='center'; tr.appendChild(td); });
    statusBody.appendChild(tr);
  }
}

/* ===== Exports ===== */
btnCSV.addEventListener('click', function(){
  if(!lastSweep){ alert('Please click Solve first.'); return; }
  const rows = [];
  rows.push(['N','BC','Re','Nu','h_Wm2K','dP_kPa','Tmax_C']);
  for(let i=0;i<lastSweep.Ns.length;i++){
    rows.push([
      lastSweep.Ns[i],
      lastSweep.bcUsed[i],
      lastSweep.Re[i],
      lastSweep.Nu[i],
      lastSweep.h[i],
      lastSweep.dPkPa[i],
      lastSweep.TmaxC[i]
    ]);
  }
  const csv = rows.map(function(r){return r.join(',');}).join('\n');
  downloadText('hx_sweep_'+Date.now()+'.csv', 'text/csv', csv);
});

btnJSON.addEventListener('click', function(){
  if(!lastSweep){ alert('Please click Solve first.'); return; }
  var inputs = {
    plate: {
      width_mm: parseFloat(plateWidth.value) || 200,
      length_mm: parseFloat(plateLength.value) || 200,
      material: plateMaterial.value,
      thickness_mm: parseFloat(plateThk.value) || 2.0
    },
    channels: {
      depth_b_mm: parseFloat(chanB.value) || 1.5,
      wall_t_mm: parseFloat(minWall.value) || 1.5,
      laminar_BC: lamBC.value
    },
    coolant: {
      type: coolantSel.value,
      flow_Lpm: parseFloat(flowLpm.value) || 0,
      inlet_C: parseFloat(inletT.value) || 25
    },
    heat_sources: heatSources
  };
  var results = [];
  for (var i = 0; i < lastSweep.Ns.length; i++){
    var N = lastSweep.Ns[i];
    var fld = (lastSweep.Tfields && lastSweep.Tfields[N]) ? lastSweep.Tfields[N] : null;
    var TminC = (fld && typeof fld.Tmin === 'number') ? fld.Tmin : null;
    results.push({
      N: N,
      BC: lastSweep.bcUsed[i],
      Re: lastSweep.Re[i],
      Nu: lastSweep.Nu[i],
      h_Wm2K: lastSweep.h[i],
      dP_kPa: lastSweep.dPkPa[i],
      Tmax_C: lastSweep.TmaxC[i],
      Tmin_C: TminC
    });
  }
  var payload = {
    exported_at: new Date().toISOString(),
    units: { h: 'W/m^2-K', dP: 'kPa', T: 'degC' },
    inputs: inputs,
    results: results
  };
  downloadText('hx_sweep_' + Date.now() + '.json', 'application/json', JSON.stringify(payload, null, 2));
});

function downloadText(filename, mime, text){
  const blob = new Blob([text], {type: mime});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

/* ===== Contour panel (heatmap rendering) ===== */
function renderContourPanel(field){
  // frame + axes
  const WmmLoc = parseNum(plateWidth,200);
  const LmmLoc = parseNum(plateLength,200);
  const panel=contourCanvas.parentElement;
  const pad=parseFloat(getComputedStyle(panel).paddingLeft)+parseFloat(getComputedStyle(panel).paddingRight);
  const innerW=panel.clientWidth - pad;

  const m={left:36,right:8,top:10,bottom:22};
  const plotWpx=Math.max(10, innerW - (m.left+m.right));
  const plotHpx=Math.max(10, plotWpx * (WmmLoc/LmmLoc));
  contourCanvas.width=Math.round(m.left+plotWpx+m.right);
  contourCanvas.height=Math.round(m.top+plotHpx+m.bottom);

  const ctx=contourCanvas.getContext('2d');
  ctx.clearRect(0,0,contourCanvas.width,contourCanvas.height);
  ctx.fillStyle='#fff'; ctx.fillRect(0,0,contourCanvas.width,contourCanvas.height);
  ctx.strokeStyle='#d9d9d9'; ctx.strokeRect(0.5,0.5,contourCanvas.width-1,contourCanvas.height-1);

  const x0c=m.left, y0c=m.top;
  const Xpx = function(x){ return x0c + (x / LmmLoc) * plotWpx; };
  const Ypx = function(y){ return (y0c + plotHpx) - (y / WmmLoc) * plotHpx; };

  // light grid
  const majorX=niceStep(LmmLoc,9), majorY=niceStep(WmmLoc,9), minorDiv=5;
  ctx.save(); ctx.lineWidth=1; ctx.strokeStyle='#eeeeee'; ctx.beginPath();
  for(let x=0;x<=LmmLoc+1e-6;x+=majorX/minorDiv){const px=Math.round(Xpx(x))+0.5;ctx.moveTo(px,y0c);ctx.lineTo(px,y0c+plotHpx);}
  for(let y=0;y<=WmmLoc+1e-6;y+=majorY/minorDiv){const py=Math.round(Ypx(y))+0.5;ctx.moveTo(x0c,py);ctx.lineTo(x0c+plotWpx,py);} ctx.stroke(); ctx.restore();
  ctx.save(); ctx.lineWidth=1; ctx.strokeStyle='#cfd4dc'; ctx.beginPath();
  for(let x=0;x<=LmmLoc+1e-6;x+=majorX){const px=Math.round(Xpx(x))+0.5;ctx.moveTo(px,y0c);ctx.lineTo(px,y0c+plotHpx);}
  for(let y=0;y<=WmmLoc+1e-6;y+=majorY){const py=Math.round(Ypx(y))+0.5;ctx.moveTo(x0c,py);ctx.lineTo(x0c+plotWpx,py);} ctx.stroke(); ctx.restore();

  // heatmap (blue->red) from field.T (ny x nx)
  const img = ctx.createImageData(Math.round(plotWpx), Math.round(plotHpx));
  const nx   = field.nx;
  const ny   = field.ny;
  const L    = field.L;
  const W    = field.W;
  const T    = field.T;
  const Tmin = field.Tmin;
  const Tmax = field.Tmax;
  const N    = field.N;
  const bc   = field.bc;

  for(let py=0; py<img.height; py++){
    const y_m = (1 - (py/(img.height-1))) * W; // bottom=0 → top=W
    const v = (y_m / W)*ny - 0.5; const j0 = Math.floor(v); const tj = clamp(v - j0, 0, 1);
    const jA = clamp(j0, 0, ny-1), jB = clamp(j0+1, 0, ny-1);
    for(let px=0; px<img.width; px++){
      const x_m = (px/(img.width-1)) * L;
      const u = (x_m / L)*nx - 0.5; const i0 = Math.floor(u); const ti = clamp(u - i0, 0, 1);
      const iA = clamp(i0, 0, nx-1), iB = clamp(i0+1, 0, nx-1);
      const T00=T[jA][iA], T10=T[jA][iB], T01=T[jB][iA], T11=T[jB][iB];
      const Ttop = T00*(1-ti) + T10*ti;
      const Tbot = T01*(1-ti) + T11*ti;
      const Tpx  = Ttop*(1-tj) + Tbot*tj;

      const s = (Tpx - Tmin) / Math.max(1e-12, (Tmax - Tmin));
      const col = jet(s);
      const r = col[0], g = col[1], b_ = col[2];
      const k = 4*(py*img.width + px);
      img.data[k+0]=r; img.data[k+1]=g; img.data[k+2]=b_; img.data[k+3]=255;
    }
  }
  ctx.putImageData(img, x0c, y0c);

  // bounding box
  ctx.strokeStyle='#0b73c8'; ctx.lineWidth=2;
  ctx.strokeRect(Math.round(x0c)+0.5, Math.round(y0c)+0.5, Math.round(plotWpx), Math.round(plotHpx));

  // axes labels
  ctx.save(); ctx.fillStyle='#111'; ctx.font='11px system-ui,sans-serif'; ctx.textAlign='center'; ctx.textBaseline='top';
  for(let x=0;x<=LmmLoc+1e-6;x+=majorX){ctx.fillText(x.toFixed(0), Xpx(x), y0c+plotHpx+4);}
  ctx.textAlign='right'; ctx.textBaseline='middle';
  for(let y=0;y<=WmmLoc+1e-6;y+=majorY){ctx.fillText(y.toFixed(0), x0c-6, Ypx(y));}
  ctx.fillStyle='#666'; ctx.textAlign='center'; ctx.fillText('x (mm)', Xpx(LmmLoc/2), y0c+plotHpx+16);
  ctx.save(); ctx.translate(x0c-26, y0c+plotHpx/2); ctx.rotate(-Math.PI/2); ctx.fillText('y (mm)', 0, 0); ctx.restore(); ctx.restore();

  // summary text
  contourSummary.textContent = 'N='+N+', BC='+bc+' • T: '+Tmin.toFixed(2)+'–'+Tmax.toFixed(2)+' °C';
}

function jet(s){
  s = clamp(s,0,1);
  // Simple JET-like (blue→cyan→yellow→red)
  const r = Math.round(255 * clamp(1.5 - Math.abs(4*s - 3), 0, 1));
  const g = Math.round(255 * clamp(1.5 - Math.abs(4*s - 2), 0, 1));
  const b = Math.round(255 * clamp(1.5 - Math.abs(4*s - 1), 0, 1));
  return [r,g,b];
}

/* ===== Events & init ===== */
function onDimsChange(){ updateChannelCounts(); drawGrid(); drawEndView(); }
['input','change'].forEach(function(evt){
  plateWidth.addEventListener(evt, onDimsChange);
  plateLength.addEventListener(evt, drawGrid);
  chanB.addEventListener(evt, drawEndView);
  minWall.addEventListener(evt, drawEndView);
  plateThk.addEventListener(evt, drawEndView);
  lamBC.addEventListener(evt, function(){ if(lastSweep) runSolve(); });
  plateMaterial.addEventListener(evt, function(){});
  coolantSel.addEventListener(evt, function(){});
  flowLpm.addEventListener(evt, function(){});
  inletT.addEventListener(evt, function(){});
});
window.addEventListener('resize', function(){ drawGrid(); drawEndView(); });

document.addEventListener('DOMContentLoaded', function(){
  updateChannelCounts(); drawGrid(); drawEndView(); refreshHSList(); updateTotalHeat();
});
