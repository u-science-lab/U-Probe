"""Self-contained, protocol-driven English HTML design reports."""
import json
import base64
from datetime import datetime
from pathlib import Path

import pandas as pd

from uprobe.core.process import filter_table
from .formatting import prepare_result_table
from .excel import column_styles, workbook_bytes


def build_report_data(df, protocol, raw_df=None, csv_filename=None):
    display_df = prepare_result_table(df)
    display_raw = prepare_result_table(raw_df) if raw_df is not None else None
    target_col = next((c for c in ("target", "gene") if c in df.columns), None)
    raw_target = next((c for c in ("target", "gene") if raw_df is not None and c in raw_df.columns), None)
    requested = []
    for target in protocol.get("targets", []):
        requested.extend(str(k) for k in target) if isinstance(target, dict) else requested.append(str(target))
    if not requested and target_col:
        requested = df[target_col].dropna().astype(str).unique().tolist()
    requested = list(dict.fromkeys(requested))
    counts = df[target_col].astype(str).value_counts().to_dict() if target_col else {}
    raw_counts = raw_df[raw_target].astype(str).value_counts().to_dict() if raw_target else {}
    filters = protocol.get("post_process", {}).get("filters", {}) or {}
    filtered = filter_table(raw_df.copy(), filters) if raw_df is not None else None
    diagnostics = []
    for key, config in filters.items():
        if not isinstance(config, dict) or "condition" not in config:
            continue
        passed = len(filter_table(raw_df.copy(), {key: config})) if raw_df is not None else None
        diagnostics.append({"metric": key, "condition": config["condition"], "passed": passed,
                            "rejected": len(raw_df) - passed if passed is not None else None})
    summary = protocol.get("summary", {}) or {}
    settings = summary.get("report", {}) or {}
    attributes = protocol.get("attributes", {}) or {}
    configured = summary.get("attributes", []) or []
    numeric = [c for c in (list(configured) or list(attributes)) if c in df and pd.api.types.is_numeric_dtype(df[c])]
    if not numeric:
        numeric = [c for c in df.select_dtypes(include="number").columns if c not in {"start", "end", "n_trans"}]
    metrics = [c for c in settings.get("key_metrics", numeric[:3]) if c in df][:3]
    table_columns = list(display_df.columns)
    statistics = []
    for col in numeric:
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        values = values[values.map(lambda x: float("-inf") < x < float("inf"))]
        if values.empty:
            continue
        unit = "°C" if attributes.get(col, {}).get("type") in {"annealing_temperature", "tm"} else ""
        histogram = pd.cut(values, bins=12, include_lowest=True).value_counts(sort=False)
        statistics.append({"name": col, "median": float(values.median()), "min": float(values.min()),
                           "max": float(values.max()), "unit": unit,
                           "bins": [{"label": str(interval), "count": int(count)} for interval, count in histogram.items()]})
    covered = sum(counts.get(t, 0) > 0 for t in requested) if target_col else None
    state = "empty" if df.empty else "partial" if covered is not None and covered < len(requested) else "complete"
    return {
        "name": str(protocol.get("name", "Probe design")), "genome": str(protocol.get("genome", "Not specified")),
        "source": str(protocol.get("extracts", {}).get("target_region", {}).get("source", "Not specified")),
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"), "state": state,
        "finalCount": len(df), "rawCount": len(raw_df) if raw_df is not None else None,
        "filterCount": len(filtered) if filtered is not None and filters else None,
        "filtersConfigured": bool(filters), "covered": covered, "targetCount": len(requested),
        "targets": [{"target": t, "raw": raw_counts.get(t, 0) if raw_target else None,
                     "final": counts.get(t, 0) if target_col else None} for t in requested],
        "diagnostics": diagnostics, "metrics": metrics, "statistics": statistics,
        "columns": list(display_df.columns), "tableColumns": table_columns,
        "columnStyles": column_styles(display_df.columns, protocol),
        "sequenceColumns": [c for c in ["target_region"] + list(protocol.get("probes", {})) if c in display_df],
        "rows": json.loads(display_df.to_json(orient="records")),
        "rawColumns": list(display_raw.columns) if display_raw is not None else [],
        "rawRows": json.loads(display_raw.to_json(orient="records")) if display_raw is not None else [],
        "csvFilename": Path(csv_filename or "probes.csv").name,
        "protocol": protocol,
    }


def save_compact_report(df, protocol, output_path, raw_df=None, csv_filename=None):
    data = build_report_data(df, protocol, raw_df, csv_filename)
    data["xlsxFilename"] = Path(data["csvFilename"]).with_suffix(".xlsx").name
    data["xlsxData"] = base64.b64encode(workbook_bytes(df, protocol)).decode("ascii")
    data["rawXlsxData"] = base64.b64encode(workbook_bytes(raw_df, protocol)).decode("ascii") if raw_df is not None else None
    # Escape HTML-sensitive characters even inside the non-executable JSON block.
    payload = json.dumps(data, ensure_ascii=True, default=str, allow_nan=False).replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(HTML.replace("$payload", payload), encoding="utf-8")
    return path


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>U-Probe Design Report</title>
<style>
:root{--ink:#25313c;--muted:#73808d;--line:#e4eaf0;--accent:#176b65}*{box-sizing:border-box}body{margin:0;background:#f7f9fb;color:var(--ink);font:14px/1.6 system-ui,-apple-system,"Segoe UI",sans-serif}main{max-width:1120px;margin:40px auto;padding:0 28px}h1{font-size:28px;letter-spacing:-.7px;margin:8px 0;font-weight:650}h2{font-size:17px;margin:0}header,.head,.actions{display:flex;align-items:center;justify-content:space-between;gap:14px}.brand{color:var(--accent);letter-spacing:.06em;font-weight:750}.muted{color:var(--muted);font-size:12px}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:16px;margin:28px 0}.card,.panel,details.section{border:1px solid var(--line);border-radius:10px;background:white}.card{padding:20px 24px}.number{font-size:36px;font-weight:650;letter-spacing:-1px;margin:3px 0}.panel{padding:24px;margin-bottom:20px}.head{margin-bottom:18px}button,input,select{font:inherit;border:1px solid var(--line);background:white;border-radius:7px;padding:8px 12px;color:var(--ink)}button,select{cursor:pointer}button:hover{border-color:var(--accent)}button:disabled{opacity:.5;cursor:default}.primary{background:var(--accent);border-color:var(--accent);color:white}.flow{display:grid;grid-template-columns:repeat(3,1fr);gap:30px}.stage{position:relative}.stage:not(:last-child):after{content:'→';position:absolute;right:-22px;top:20px;color:#a4afbb}.count{font-size:26px;font-weight:600}.track{height:5px;background:#edf1f4;border-radius:4px;margin-top:12px}.bar{height:5px;background:#83b4ad;border-radius:4px}.quality{display:flex;gap:32px;flex-wrap:wrap;margin-top:22px;border-top:1px solid var(--line);padding-top:18px}.quality b{display:block;font-weight:550}.scroll{overflow:auto;max-width:100%}table{width:100%;border-collapse:collapse;font-size:13px;text-align:left}td,th{padding:12px 14px;border-bottom:1px solid var(--line);white-space:nowrap}th{color:var(--muted);font-size:12px;font-weight:550;background:#f8fafb;cursor:pointer}tbody tr:last-child td{border-bottom:0}.mono-cell{font-family:ui-monospace,Consolas,monospace}.sequence{background:#edf5f3;box-shadow:inset 1px 0 #d6e5e1,inset -1px 0 #d6e5e1;font-family:ui-monospace,Consolas,monospace}.sequence.target-sequence{background:#eef2fa;box-shadow:inset 1px 0 #dce3f0,inset -1px 0 #dce3f0}.good{color:var(--accent)}.warning{color:#99661a}.notice{border:1px solid #eddbb8;background:#fff9ee;padding:18px 22px;border-radius:9px;margin:24px 0;color:#81591e}.notice button{margin-top:10px}.hidden{display:none!important}details.section{padding:18px 24px;margin-bottom:12px}summary{cursor:pointer;font-weight:550}summary .muted{float:right;font-weight:400}.body{margin-top:20px}.pager{display:flex;justify-content:space-between;align-items:center;margin-top:16px}.pager button{margin-left:8px}.mono{font:12px/1.6 ui-monospace,Consolas,monospace;overflow-wrap:anywhere;white-space:pre-wrap;max-width:800px}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;background:#f7f9fb;padding:16px;border-radius:7px}.chart{display:flex;align-items:flex-end;gap:6px;height:150px;border-bottom:1px solid var(--line);margin-top:20px}.chart i{flex:1;background:#83b4ad;border-radius:3px 3px 0 0;min-height:1px}.axis{display:flex;justify-content:space-between;color:var(--muted);font-size:12px;margin-top:8px}footer{color:var(--muted);font-size:12px;padding:22px 0}dialog{border:1px solid var(--line);border-radius:12px;width:min(980px,96vw);padding:24px;color:var(--ink)}dialog::backdrop{background:#17323b66}.close{float:right}@media(max-width:700px){main{padding:0 16px;margin:20px auto}header{align-items:flex-start;flex-direction:column}.cards{gap:8px}.card{padding:12px}.number{font-size:28px}.panel{padding:18px}.quality{gap:16px}.head{flex-wrap:wrap}summary .muted{float:none;display:block}}@media print{body{background:white}button,input,select,.pager{display:none}main{max-width:none;margin:0}.panel,.card{break-inside:avoid}}
</style></head><body><main>
<header><div><div class="brand">U-PROBE / DESIGN REPORT</div><h1 id="name"></h1><div class="muted" id="subtitle"></div></div><div class="actions"><button id="viewRaw" class="hidden">View raw</button><button class="primary" id="download">Download XLSX</button><button id="print">Print</button></div></header>
<div class="notice hidden" id="notice"><strong id="noticeTitle"></strong><div id="noticeText"></div></div>
<div class="cards"><div class="card"><div class="muted">Final probes</div><div class="number" id="final"></div><div class="muted" id="finalNote"></div></div><div class="card"><div class="muted">Target coverage</div><div class="number" id="coverage"></div><div class="muted" id="coverageNote"></div></div><div class="card"><div class="muted">Filter pass rate</div><div class="number" id="rate"></div><div class="muted" id="rateNote"></div></div></div>
<section class="panel"><div class="head"><h2>Selection overview</h2><span class="muted">Candidates to final results</span></div><div class="flow" id="flow"></div><div class="quality" id="quality"></div></section>
<section class="panel"><div class="head"><h2>Target summary</h2><span class="muted">Targets with at least one final probe</span></div><div class="scroll" id="targets"></div></section>
<section class="panel"><div class="head"><div><h2>Final probes</h2><div class="muted">Scroll horizontally to view all fields; click a column to sort. Matching colors group sequences and their attributes.</div></div><input id="search" aria-label="Search final probes" placeholder="Search probes or targets"></div><div class="scroll" id="results"></div><div class="pager"><span class="muted" id="pageInfo"></span><div><button id="prev">Previous</button><button id="next">Next</button></div></div></section>
<details class="section" id="diagnostics"><summary>Filter diagnostics <span class="muted">Conditions and independent pass counts</span></summary><div class="body"><p class="muted">Each condition is checked independently on raw candidates. Rejection counts can overlap and must not be added together. Conditions for missing columns are skipped by the pipeline.</p><div class="scroll" id="filterTable"></div></div></details>
<details class="section"><summary>Quality statistics <span class="muted">One metric at a time</span></summary><div class="body"><select id="metric" aria-label="Quality metric"></select><p class="muted" id="statSummary"></p><div class="chart" id="chart"></div><div class="axis"><span id="low"></span><span id="high"></span></div><p class="muted">Final probes only. Hover over a bar for its interval and count.</p></div></details>
<footer>U-Probe · Final results and unfiltered candidates are reported separately. Identifiers and protocol values are retained as supplied.</footer>
</main><dialog id="rawDialog"><button class="close" id="closeRaw">Close</button><h2>Raw candidates</h2><p class="notice">Unfiltered candidates for troubleshooting. They are not passing probes. Preview shows up to 100 rows.</p><button id="downloadRaw">Download raw XLSX</button><div class="scroll" id="rawTable" style="max-height:55vh;margin-top:16px"></div></dialog>
<script type="application/json" id="report-data">$payload</script>
<script>
const data=JSON.parse(document.getElementById('report-data').textContent);
const el=id=>document.getElementById(id),text=(id,value)=>el(id).textContent=value;
const esc=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const cellClass=c=>{const s=data.columnStyles[c];return s?' style="background:#'+s.background+'"'+(s.sequence?' class="mono-cell"':''):''};
const label=value=>String(value).replace(/[_.]/g,' ').replace(/\b\w/g,c=>c.toUpperCase());
const num=value=>value===null||value===undefined?'—':(Math.abs(Number(value))<0.005?0:Number(value)).toLocaleString('en-US',{maximumFractionDigits:2});
function table(columns,rows){return '<table><thead><tr>'+columns.map(c=>'<th>'+esc(label(c))+'</th>').join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+columns.map(c=>'<td'+cellClass(c)+'>'+esc(typeof r[c]==='object'&&r[c]!==null?JSON.stringify(r[c]):r[c]??'—')+'</td>').join('')+'</tr>').join('')+'</tbody></table>'}
function downloadWorkbook(encoded,filename){const bytes=Uint8Array.from(atob(encoded),c=>c.charCodeAt(0));const url=URL.createObjectURL(new Blob([bytes],{type:'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'}));const a=document.createElement('a');a.href=url;a.download=filename;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}
text('name',data.name);document.title='U-Probe · '+data.name;text('subtitle',data.genome+' · '+data.source+' · '+data.generated);
text('final',num(data.finalCount));text('finalNote',data.filtersConfigured?'After filtering and post-processing':'No filter conditions configured');
text('coverage',data.targetCount&&data.covered!==null?num(data.covered/data.targetCount*100)+'%':'—');text('coverageNote',data.covered!==null?num(data.covered)+' / '+num(data.targetCount)+' requested targets':'Target information unavailable');
text('rate',data.filtersConfigured&&data.rawCount>0&&data.filterCount!==null?num(data.filterCount/data.rawCount*100)+'%':'—');text('rateNote',!data.filtersConfigured?'Filtering not configured':data.filterCount===null?'Raw candidate counts unavailable':num(data.filterCount)+' / '+num(data.rawCount)+' candidates passed all conditions');
const stages=[['Raw candidates',data.rawCount],['Condition filters',data.filtersConfigured?data.filterCount:null],['Final probes',data.finalCount]];el('flow').innerHTML=stages.map(([name,count])=>'<div class="stage"><div class="muted">'+name+'</div><div class="count">'+num(count)+'</div><div class="track"><div class="bar" style="width:'+(data.rawCount>0&&count!==null?Math.min(100,count/data.rawCount*100):0)+'%"></div></div></div>').join('');
el('quality').innerHTML=data.metrics.map(m=>data.statistics.find(s=>s.name===m)).filter(Boolean).map(s=>'<div><span class="muted">'+esc(label(s.name))+'</span><b>'+num(s.median)+' '+esc(s.unit)+'</b><span class="muted">Median · '+num(s.min)+'–'+num(s.max)+'</span></div>').join('')||'<span class="muted">No numeric quality metrics available.</span>';
el('targets').innerHTML=data.targets.length?table(['target','raw_candidates','final_probes','status'],data.targets.map(t=>({target:t.target,raw_candidates:t.raw,final_probes:t.final,status:t.final===null?'Unavailable':t.final>0?'Covered':'No result'}))):'<p class="muted">Target information unavailable.</p>';
if(data.state!=='complete'){el('notice').classList.remove('hidden');text('noticeTitle',data.state==='empty'?'No final probes':'Some targets have no final probes');text('noticeText',data.state==='empty'?'No candidates remain after selection. Review filter diagnostics and raw candidates.':'Review targets with no result before using this design.');el('diagnostics').open=data.state==='empty'}
if(data.rawRows.length)el('viewRaw').classList.remove('hidden');
el('viewRaw').onclick=()=>{el('rawTable').innerHTML=table(data.rawColumns,data.rawRows.slice(0,100));el('rawDialog').showModal()};el('closeRaw').onclick=()=>el('rawDialog').close();el('downloadRaw').onclick=()=>downloadWorkbook(data.rawXlsxData,data.xlsxFilename.replace(/\.xlsx$/,'_raw.xlsx'));
el('download').disabled=!data.rows.length;el('download').onclick=()=>downloadWorkbook(data.xlsxData,data.xlsxFilename);el('print').onclick=()=>window.print();
el('filterTable').innerHTML=data.diagnostics.length?table(['metric','condition','passed','rejected'],data.diagnostics):'<p class="muted">No filter conditions configured.</p>';
let page=0,sortColumn=null,ascending=true;function renderRows(){const query=el('search').value.toLowerCase();let rows=data.rows.filter(r=>data.columns.some(c=>String(r[c]??'').toLowerCase().includes(query)));if(sortColumn)rows.sort((a,b)=>{const av=a[sortColumn],bv=b[sortColumn];return (typeof av==='number'&&typeof bv==='number'?av-bv:String(av??'').localeCompare(String(bv??'')))*(ascending?1:-1)});const pages=Math.max(1,Math.ceil(rows.length/20));page=Math.min(page,pages-1);const visible=rows.slice(page*20,page*20+20);el('results').innerHTML='<table><thead><tr>'+data.tableColumns.map(c=>'<th data-column="'+esc(c)+'">'+esc(label(c))+'</th>').join('')+'</tr></thead><tbody>'+visible.map(r=>'<tr>'+data.tableColumns.map(c=>'<td'+cellClass(c)+'>'+esc(typeof r[c]==='object'&&r[c]!==null?JSON.stringify(r[c]):r[c]??'—')+'</td>').join('')+'</tr>').join('')+'</tbody></table>'+(rows.length?'':'<p class="muted">No probes to display.</p>');text('pageInfo',num(rows.length)+' probes · Page '+(page+1)+' / '+pages);el('prev').disabled=page===0;el('next').disabled=page>=pages-1;el('results').querySelectorAll('[data-column]').forEach(th=>th.onclick=()=>{const c=th.dataset.column;ascending=sortColumn===c?!ascending:true;sortColumn=c;page=0;renderRows()})}
el('search').oninput=()=>{page=0;renderRows()};el('prev').onclick=()=>{page--;renderRows()};el('next').onclick=()=>{page++;renderRows()};renderRows();
el('metric').innerHTML=data.statistics.map((s,i)=>'<option value="'+i+'">'+esc(label(s.name))+'</option>').join('');function histogram(){const s=data.statistics[Number(el('metric').value)];if(!s){text('statSummary','No statistics available.');return}text('statSummary','Median '+num(s.median)+' '+s.unit+' · Range '+num(s.min)+'–'+num(s.max));const max=Math.max(...s.bins.map(b=>b.count),1);el('chart').innerHTML=s.bins.map(b=>'<i title="'+esc(b.label+': '+b.count+' probes')+'" style="height:'+b.count/max*140+'px"></i>').join('');text('low',num(s.min)+' '+s.unit);text('high',num(s.max)+' '+s.unit)}el('metric').onchange=histogram;histogram();
</script></body></html>'''
