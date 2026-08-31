import React, { useMemo, useState } from 'react'
import { scoreBoth, differential, getApiBase } from './api/client'

const API_BASE = getApiBase()

// Small default fixtures for demo (mirrors tests: injected signal)
const DEFAULT_GENE_SETS: Record<string, string[]> = {
  HALLMARK_INTERFERON_RESPONSE: ['IFIT1', 'MX1', 'ISG15', 'OAS1', 'RSAD2', 'IRF7'],
  HALLMARK_TNFA_VIA_NFKB: ['NFKB1', 'TNFAIP3', 'IL6', 'CXCL2', 'BCL2A1'],
  HALLMARK_P53_PATHWAY: ['CDKN1A', 'MDM2', 'BAX', 'GADD45A', 'TP53'],
}

const DEFAULT_EXPRESSION_CSV = `sample,gene,group
S1,IFIT1,10
S1,MX1,9.5
S1,ISG15,9
S1,OAS1,8.5
S1,RSAD2,9.2
S1,IRF7,8.8
S1,NFKB1,4
S1,TNFAIP3,3.5
S1,CDKN1A,5
S2,IFIT1,9.8
S2,MX1,9.2
S2,ISG15,9.1
S2,OAS1,8.7
S2,RSAD2,9.0
S2,IRF7,8.9
S2,NFKB1,4.2
S2,CDKN1A,5.1
S3,IFIT1,3
S3,MX1,3.2
S3,ISG15,2.9
S3,OAS1,3.1
S3,RSAD2,2.8
S3,IRF7,3.0
S3,NFKB1,4.5
S3,TNFAIP3,4.0
S3,CDKN1A,5.2
S4,IFIT1,2.8
S4,MX1,2.9
S4,ISG15,3.0
S4,OAS1,2.7
S4,RSAD2,3.1
S4,IRF7,2.9
S4,NFKB1,4.6
S4,CDKN1A,5.0`

function parseExpressionCSV(csv: string): { expr: Record<string, Record<string, number>>, groups: Record<string, string> } {
  const lines = csv.trim().split('\n')
  const header = lines[0].split(',').map(s => s.trim())
  const sIdx = header.indexOf('sample')
  const gIdx = header.indexOf('gene')
  const grpIdx = header.indexOf('group')
  // Support both long format (sample,gene,value,group) and wide? Here long
  // Expect columns: sample, gene, value OR sample,gene,group where value is next? Actually default uses sample,gene,group? No value column - we treat third as value and look for group in header? S1,IFIT1,10
  // Detect value column: if header has 'value' else assume third column is value
  let vIdx = header.indexOf('value')
  if (vIdx === -1) {
    // fallback: if we have 3 cols and second col is gene, third is numeric value (as in default truncated)
    // Let's detect type of third column
    vIdx = 2
  }
  const expr: Record<string, Record<string, number>> = {}
  const groups: Record<string, string> = {}
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(',').map(s => s.trim())
    if (parts.length < 3) continue
    const sample = parts[sIdx]
    const gene = parts[gIdx]
    const val = parseFloat(parts[vIdx])
    if (isNaN(val)) continue
    if (!expr[sample]) expr[sample] = {}
    expr[sample][gene] = val
    if (grpIdx !== -1 && parts[grpIdx] && grpIdx !== vIdx) {
      groups[sample] = parts[grpIdx]
    }
  }
  // If groups not parsed, infer from sample name pattern: S1,S2 are group A etc for demo
  if (Object.keys(groups).length === 0) {
    for (const s of Object.keys(expr)) {
      // S1,S2 high interferon -> group A, S3,S4 low -> group B
      groups[s] = (s === 'S1' || s === 'S2') ? 'A' : 'B'
    }
  }
  return { expr, groups }
}

function parseWideCSV(csv: string): { expr: Record<string, Record<string, number>>, groups: Record<string, string> } {
  // Wide format: first row header gene names, first column sample
  // Optional second file for groups, but we try to parse group suffix: sample_A etc
  const lines = csv.trim().split('\n').filter(Boolean)
  if (lines.length < 2) return { expr: {}, groups: {} }
  const header = lines[0].split(',').map(s => s.trim())
  // header[0] is sample, rest are genes
  const genes = header.slice(1)
  const expr: Record<string, Record<string, number>> = {}
  const groups: Record<string, string> = {}
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split(',').map(s => s.trim())
    const sample = parts[0]
    if (!sample) continue
    expr[sample] = {}
    for (let j = 0; j < genes.length; j++) {
      const v = parseFloat(parts[j + 1])
      if (!isNaN(v)) expr[sample][genes[j]] = v
    }
    // group hint: sample contains _A or _B or look up not provided
    if (sample.includes('_A')) groups[sample] = 'A'
    else if (sample.includes('_B')) groups[sample] = 'B'
  }
  return { expr, groups }
}

export default function App() {
  const [geneSetsText, setGeneSetsText] = useState(JSON.stringify(DEFAULT_GENE_SETS, null, 2))
  const [exprText, setExprText] = useState(DEFAULT_EXPRESSION_CSV)
  const [format, setFormat] = useState<'long' | 'wide'>('long')
  const [groupsText, setGroupsText] = useState('S1,A\nS2,A\nS3,B\nS4,B')
  const [loading, setLoading] = useState(false)
  const [result, setResult] = useState<any>(null)
  const [diff, setDiff] = useState<any>(null)
  const [error, setError] = useState<string | null>(null)
  const [health, setHealth] = useState<any>(null)

  const parsedGeneSets = useMemo(() => {
    try { return JSON.parse(geneSetsText) as Record<string, string[]> } catch { return null }
  }, [geneSetsText])

  const run = async () => {
    setError(null); setLoading(true); setResult(null); setDiff(null)
    try {
      if (!parsedGeneSets) throw new Error('Gene sets JSON invalid')
      let expr: Record<string, Record<string, number>> = {}
      let groups: Record<string, string> = {}
      if (format === 'long') {
        // Try long parser; if header contains gene count wide, fallback to wide
        const header = exprText.trim().split('\n')[0]
        if (header.split(',').length > 3 && !header.includes('value')) {
          const w = parseWideCSV(exprText)
          expr = w.expr
          groups = w.groups
        } else {
          // Our default CSV is actually sample,gene,expression where third column numeric, no header value -> parse manually
          // Custom quick long with numeric third column
          const lines = exprText.trim().split('\n')
          const h = lines[0].split(',').map(s=>s.trim().toLowerCase())
          if (h[0]==='sample' && h[1]==='gene') {
            // long
            // If header is sample,gene,group but rows are sample,gene,value, the third header is misleading. Detect numeric in row 1.
            const firstRow = lines[1]?.split(',')
            const thirdIsNum = firstRow ? !isNaN(parseFloat(firstRow[2])) : false
            if (thirdIsNum) {
              // treat third column as value, ignore group header
              const tmpExpr: Record<string, Record<string, number>> = {}
              for (let i=1;i<lines.length;i++) {
                const p = lines[i].split(',').map(s=>s.trim())
                if (p.length<3) continue
                const [s,g,v] = p
                const num = parseFloat(v)
                if (isNaN(num)) continue
                if (!tmpExpr[s]) tmpExpr[s]={}
                tmpExpr[s][g]=num
              }
              expr = tmpExpr
              // groups from separate textarea
              groups = Object.fromEntries(groupsText.trim().split('\n').filter(Boolean).map(l=>{const [k,v]=l.split(',').map(s=>s.trim()); return [k,v]}))
            } else {
              const p = parseExpressionCSV(exprText)
              expr = p.expr
              groups = p.groups
            }
          } else {
            const p = parseWideCSV(exprText)
            expr = p.expr
            groups = p.groups
          }
          // if groupsText overrides, use it
          if (groupsText.trim()) {
            const g2 = Object.fromEntries(groupsText.trim().split('\n').filter(Boolean).map(l=>{const [k,v]=l.split(',').map(s=>s.trim()); return [k,v]}))
            if (Object.keys(g2).length) groups = { ...groups, ...g2 }
          }
        }
      } else {
        const w = parseWideCSV(exprText)
        expr = w.expr
        if (groupsText.trim()) {
          groups = Object.fromEntries(groupsText.trim().split('\n').filter(Boolean).map(l=>{const [k,v]=l.split(',').map(s=>s.trim()); return [k,v]}))
        } else {
          groups = w.groups
        }
      }
      if (!Object.keys(expr).length) throw new Error('No expression data parsed')
      const res = await scoreBoth(expr, parsedGeneSets)
      setResult(res)
      // differential using ssGSEA scores
      const diffReq = { scores: res.ssgsea, groups }
      // if groups has at least 2 unique labels
      const uniq = [...new Set(Object.values(groups))]
      if (uniq.length >= 2) {
        const d = await differential(diffReq.scores, diffReq.groups)
        setDiff(d)
      }
    } catch (e:any) {
      setError(e.message || String(e))
    } finally { setLoading(false) }
  }

  const checkHealth = async () => {
    try {
      const r = await fetch(`${API_BASE}/health`)
      const j = await r.json()
      setHealth(j)
    } catch {
      // fallback: try relative health if API_BASE was absolute and failed
      try { const r = await fetch('/health'); const j = await r.json(); setHealth(j) } catch {}
    }
  }

  // Compute top pathway for chart
  const topPathways = useMemo(() => {
    if (!diff) return []
    return [...diff].sort((a:any,b:any)=> a.q_value - b.q_value).slice(0, 6)
  }, [diff])

  return (
    <div>
      <div className="header">
        <h1>Pathway Activity Dashboard — ssGSEA (Barbie 2009) · Combined z-score (Lee 2008)</h1>
        <p>
          Context-specific pathway activity inference from expression. Paste a matrix (samples × genes), provide gene sets (GMT-style or JSON), compute per-sample ssGSEA and z-score activity, view ranked differential activity (Wilcoxon + BH-FDR) and method agreement (Spearman/Pearson). Scoring runs ungated; curated artifact endpoints require <span className="code">MODEL_RELEASE_APPROVED=true</span>.
        </p>
        <div className="meta">
          <span className="badge">React + Vite + TypeScript</span>
          <span className="badge">FastAPI backend</span>
          <span className="badge">MSigDB Hallmark · Reactome GMT</span>
          <span className="badge">Fail-closed release gate</span>
        </div>
      </div>

      <div className="container">
        <div className="grid">
          <div className="card">
            <h2>1 · Expression matrix</h2>
            <div className="muted">Paste CSV. Supports long (sample,gene,value) or wide (sample × gene). Demo has injected signal: IFN genes high in S1/S2.</div>
            <div style={{ display: 'flex', gap: 8, margin: '8px 0' }}>
              <select value={format} onChange={e=>setFormat(e.target.value as any)} style={{ width: 140 }}>
                <option value="long">Long format</option>
                <option value="wide">Wide (samples × genes)</option>
              </select>
              <button className="btn ghost" onClick={checkHealth}>Check /health</button>
            </div>
            <textarea rows={10} value={exprText} onChange={e=>setExprText(e.target.value)} placeholder="sample,gene,value" />
            <h3>Sample groups (sample,group per line)</h3>
            <textarea rows={4} value={groupsText} onChange={e=>setGroupsText(e.target.value)} />
            {health && <pre className="muted" style={{ background: '#f1f5f9', padding: 8, borderRadius: 8, marginTop: 8, overflow: 'auto' }}>{JSON.stringify(health, null, 2)}</pre>}
          </div>

          <div className="card">
            <h2>2 · Gene sets (GMT / JSON)</h2>
            <div className="muted">JSON: {"{ pathway: [genes] }"} • Or upload GMT file to backend via <span className="code">--gmt-path</span>. Uses MSigDB Hallmark (50) or Reactome GMT.</div>
            <textarea rows={14} value={geneSetsText} onChange={e=>setGeneSetsText(e.target.value)} />
            <div style={{ display:'flex', gap:8, marginTop:10 }}>
              <button className="btn" onClick={run} disabled={loading || !parsedGeneSets}>
                {loading ? 'Scoring…' : 'Run ssGSEA + z-score'}
              </button>
              <button className="btn secondary" onClick={()=>{setGeneSetsText(JSON.stringify(DEFAULT_GENE_SETS,null,2)); setExprText(DEFAULT_EXPRESSION_CSV); setGroupsText('S1,A\nS2,A\nS3,B\nS4,B')}}>Reset demo</button>
            </div>
            {!parsedGeneSets && <div className="alert error" style={{marginTop:8}}>Invalid JSON</div>}
            {error && <div className="alert error" style={{marginTop:8}}>{error}</div>}
            <div className="muted" style={{marginTop:8}}>Scoring endpoints are <b>ungated</b> (deterministic stats). Curated endpoints (<span className="code">/api/v1/curated/*</span>) require <span className="code">MODEL_RELEASE_APPROVED=true</span> + <span className="code">APPROVED_ARTIFACT_REVISION</span>.</div>
          </div>
        </div>

        {result && (
          <>
            <div className="card" style={{ marginTop: 16 }}>
              <h2>Correlation — ssGSEA vs z-score (per pathway)</h2>
              <div className="muted">Spearman and Pearson per pathway quantify method agreement on this cohort.</div>
              <table className="table" style={{ marginTop: 8 }}>
                <thead><tr><th>Pathway</th><th>Spearman r</th><th>Pearson r</th><th>n</th></tr></thead>
                <tbody>
                  {result.correlation?.map((c:any)=>(
                    <tr key={c.pathway}>
                      <td><b>{c.pathway}</b></td>
                      <td>{c.spearman_r?.toFixed(3)}</td>
                      <td>{c.pearson_r?.toFixed(3)}</td>
                      <td>{c.n_samples}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>

            <div className="grid" style={{ marginTop: 16 }}>
              <div className="card">
                <h2>Top differentially active pathways</h2>
                <div className="muted">Wilcoxon rank-sum per pathway (group B vs A) + BH-FDR. Bar is -log10(q).</div>
                {diff ? (
                  <>
                    <div style={{ display: 'flex', flexDirection: 'column', gap: 8, marginTop: 10 }}>
                      {topPathways.map((d:any)=>{
                        const q = d.q_value || 1
                        const neglog = q>0 ? Math.min(6, -Math.log10(q)) : 6
                        const width = (neglog/6)*100
                        return (
                          <div key={d.pathway}>
                            <div style={{ display:'flex', justifyContent:'space-between', fontSize:11 }}>
                              <b>{d.pathway}</b>
                              <span className="muted">Δ {d.delta?.toFixed(2)} · q {q.toExponential(1)}</span>
                            </div>
                            <div style={{ background:'#e2e8f0', borderRadius:999, height:14, marginTop:4, overflow:'hidden' }}>
                              <div className="bar" style={{ width:`${width}%`, height:'100%' }} />
                            </div>
                          </div>
                        )
                      })}
                    </div>
                    <table className="table" style={{ marginTop: 12 }}>
                      <thead><tr><th>Pathway</th><th>Mean A</th><th>Mean B</th><th>Δ</th><th>p</th><th>q</th></tr></thead>
                      <tbody>
                        {diff.map((d:any)=>(
                          <tr key={d.pathway}>
                            <td>{d.pathway}</td>
                            <td>{d.mean_A?.toFixed(2)}</td>
                            <td>{d.mean_B?.toFixed(2)}</td>
                            <td style={{ color: d.delta>0 ? '#059669' : '#dc2626', fontWeight:700 }}>{d.delta?.toFixed(2)}</td>
                            <td>{d.p_value?.toExponential(1)}</td>
                            <td><span className={`pill ${d.q_value < 0.05 ? 'sig' : 'ns'}`}>{d.q_value?.toExponential(1)}</span></td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </>
                ) : <div className="muted">No groups supplied or groups not distinct.</div>}
              </div>

              <div className="card">
                <h2>Per-sample pathway activity (ssGSEA)</h2>
                <div className="muted">Heatmap (blue→low, red→high). Z-score view available via toggle in API response.</div>
                {result.ssgsea && (()=> {
                  const samples = result.samples as string[]
                  const pathways = result.pathways as string[]
                  const scores = result.ssgsea as Record<string, Record<string, number>>
                  // compute min/max for color scale
                  const vals = samples.flatMap(s=> pathways.map(p=> scores[s][p])).filter(v=> typeof v==='number')
                  const min = Math.min(...vals), max = Math.max(...vals)
                  const color = (v:number)=> {
                    const t = (v - min)/ (max - min || 1)
                    // diverging: low #3b82f6 -> high #ef4444 via white mid
                    const r = Math.round(59 + t*180)
                    const g = Math.round(130 - t*60)
                    const b = Math.round(246 - t*180)
                    return `rgb(${r},${g},${b})`
                  }
                  return (
                    <>
                      <div style={{ display:'grid', gridTemplateColumns:`120px repeat(${samples.length},1fr)`, gap:4, fontSize:11, marginTop:8 }}>
                        <div></div>
                        {samples.map(s=> <div key={s} style={{ textAlign:'center', fontWeight:700, fontSize:10 }}>{s}</div>)}
                        {pathways.map(p=> (
                          <React.Fragment key={p}>
                            <div style={{ fontWeight:600, fontSize:10, overflow:'hidden', textOverflow:'ellipsis', whiteSpace:'nowrap' }}>{p}</div>
                            {samples.map(s=> {
                              const v = scores[s][p]
                              return <div key={s} className="heat-cell" style={{ background: color(v) }} title={`${v.toFixed(2)}`}>{v.toFixed(1)}</div>
                            })}
                          </React.Fragment>
                        ))}
                      </div>
                      <div className="stat-row" style={{ marginTop: 12 }}>
                        <div className="stat"><div className="label">Samples</div><div className="value">{samples.length}</div></div>
                        <div className="stat"><div className="label">Pathways</div><div className="value">{pathways.length}</div></div>
                        <div className="stat"><div className="label">Score range</div><div className="value" style={{fontSize:13}}>{min.toFixed(1)} → {max.toFixed(1)}</div></div>
                      </div>
                    </>
                  )
                })()}
              </div>
            </div>
          </>
        )}

        <div className="muted" style={{ textAlign:'center', marginTop: 20 }}>
          Citations: Barbie et al. Nature 2009 (ssGSEA) · Lee et al. PLoS Comp Biol 2008 (combined z-score) · Subramanian et al. PNAS 2005 (GSEA) · Hänzelmann et al. BMC Bioinf 2013 (GSVA). Data sources: MSigDB Hallmark (gsea-msigdb.org, 50 hallmarks, non-commercial) · Reactome (reactome.org/download-data, CC BY 4.0). Synthetic fixture below validates pipeline recovers injected signal; not a biological finding.
        </div>
      </div>
    </div>
  )
}
