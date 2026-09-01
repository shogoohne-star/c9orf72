// Chromium で SVG を PNG に焼く（ライトテーマ固定、2 倍解像度）
const fs = require('fs');
const path = require('path');
const { chromium } = require(process.env.PW || '/opt/node22/lib/node_modules/playwright');

const css = fs.readFileSync('style.css', 'utf8');
const OUT = 'png';
fs.mkdirSync(OUT, { recursive: true });

const CHARTS = ['chart_allowed','chart_strain','chart_em','chart_hit','chart_onoff',
                'chart_pos','chart_gap','chart_robust','fig_anchor','fig_nac','fig_em'];
const MOLS = ['A_tz','A_cp','B_tz','B_cp','C_tz','C_cp'];

function clean(s){
  return s.replace(/<\?xml[^>]*\?>\s*/,'').replace(/<!--[\s\S]*?-->/g,'');
}

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 3 });

  for (const name of CHARTS.concat(MOLS)) {
    let svg = clean(fs.readFileSync(`assets/${name}.svg`, 'utf8'));
    const vb = /viewBox=['"]([\d.\s-]+)['"]/.exec(svg);
    const [, , w, h] = vb[1].trim().split(/\s+/).map(Number);
    const isMol = MOLS.includes(name);
    if (isMol) svg = svg.replace(/\s(width|height)='[^']*'/g, '');
    const html = `<!doctype html><html data-theme="light"><head><meta charset="utf-8">
<style>${css}
html,body{margin:0;padding:0;background:#FFFFFF}
#box{width:${w}px;height:${h}px;background:#FFFFFF;color:#19191E}
#box > svg{width:${w}px;height:${h}px;display:block}
figure{all:unset}</style></head><body><div id="box">${svg}</div></body></html>`;
    await page.setViewportSize({ width: Math.ceil(w) + 20, height: Math.ceil(h) + 20 });
    await page.setContent(html, { waitUntil: 'load' });
    await page.evaluate(() => document.fonts.ready);
    await page.locator('#box').screenshot({ path: path.join(OUT, `${name}.png`) });
    console.log(name, `${w}x${h}`);
  }
  await browser.close();
})();
