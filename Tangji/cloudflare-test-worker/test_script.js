
// Mock DOM
global.document = {
  getElementById: function(id) {
    if (id === 'n-ca4') return { value: '0.022' };
    if (id === 'n-a2') return { value: '0.60' };
    if (id === 'n-text') return { innerHTML: '' };
    if (id === 'n-gauge') return {};
    return null;
  }
};
global.echarts = {
  init: function() {
    return {
      setOption: function(opt) {
        console.log('setOption called');
      }
    };
  }
};

      function calcNom() {
        const ca4 = parseFloat(document.getElementById('n-ca4').value) || 0;
        const a2 = parseFloat(document.getElementById('n-a2').value) || 0;
        const logit = -2.5 + 3.2 * ca4 - 2.0 * a2;
        const prob = 1 / (1 + Math.exp(-logit));
        const pct = (prob * 100).toFixed(1);
        let level = prob < 0.3 ? 'Low' : prob < 0.6 ? 'Intermediate' : 'High';
        let color = prob < 0.3 ? '#00a65a' : prob < 0.6 ? '#f39c12' : '#dd4b39';
        document.getElementById('n-text').innerHTML = '<span style="color:' + color + '">' + level + ' Risk: ' + pct + '%</span>';
        const chart = echarts.init(document.getElementById('n-gauge'));
        chart.setOption({
          series: [{
            type: 'gauge',
            min: 0, max: 100,
            axisLine: { lineStyle: { width: 20, color: [[0.3, '#00a65a'], [0.6, '#f39c12'], [1, '#dd4b39']] } },
            pointer: { length: '60%', width: 5 },
            detail: { fontSize: 20, formatter: '{value}%' },
            data: [{ value: parseFloat(pct), name: 'Risk %' }]
          }]
        });
      }
      calcNom();
    