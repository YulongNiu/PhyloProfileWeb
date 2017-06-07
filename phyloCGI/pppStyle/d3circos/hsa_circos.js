var gieStainColor = {
    gpos100: 'rgb(0,0,0)',
    gpos: 'rgb(0,0,0)',
    gpos75: 'rgb(130,130,130)',
    gpos66: 'rgb(160,160,160)',
    gpos50: 'rgb(200,200,200)',
    gpos33: 'rgb(210,210,210)',
    gpos25: 'rgb(200,200,200)',
    gvar: 'rgb(220,220,220)',
    gneg: 'rgb(255,255,255)',
    acen: 'rgb(217,47,39)',
    stalk: 'rgb(100,127,164)',
    select: 'rgb(135,177,255)'
}

var drawCircos = function (error, genome, cytobands, linkage, freqEu, freqArc, freqBac) {
    var width = 800
    var circos = new Circos({
        container: '#pppcircos',
        width: width,
        height: width
    })

    cytobands = cytobands
        .map(function (d) {
        return {
            block_id: d.chrom,
            start: parseInt(d.chromStart),
            end: parseInt(d.chromEnd),
            gieStain: d.gieStain,
            name: d.name
        }
    })

    linkage = linkage
        .map(function (d) {
            return {
                source: {
                    id: d.source_id,
                    start: parseInt(d.source_start),
                    end: parseInt(d.source_end) + 1500000,
                    position: (parseInt(d.source_start) + parseInt(d.source_end)) / 2,
                    label: d.source_label
                },
                target: {
                    id: d.target_id,
                    start: parseInt(d.target_start),
                    end: parseInt(d.target_end) + 1500000,
                    position: (parseInt(d.target_start) + parseInt(d.target_end)) / 2,
                    label: d.target_label
                },
                color: d.linkage_color,
                jaccard: d.Jaccard,
                cor: d.Cor
            }
        })

    freqEu = freqEu.map(function (d) {
        return {
            block_id: d.chromosome,
            position: (parseInt(d.start) + parseInt(d.end)) / 2,
            value: parseFloat(d.value),
            label: d.label
        }
    })

    freqArc = freqArc.map(function (d) {
        return {
            block_id: d.chromosome,
            position: (parseInt(d.start) + parseInt(d.end)) / 2,
            value: parseFloat(d.value),
            label: d.label
        }
    })

    freqBac = freqBac.map(function (d) {
        return {
            block_id: d.chromosome,
            position: (parseInt(d.start) + parseInt(d.end)) / 2,
            value: parseFloat(d.value),
            label: d.label
        }
    })

    circos
        .layout(
            genome,
            {
                innerRadius: width/2 - 150,
                outerRadius: width/2 - 130,
                gap: 0.02,
                labels: {
                    size: 8,
                    radialOffset: 140
                },
                ticks: {
                    spacing: 5000000,
                    display: true,
                    labelDenominator: 1000000,
                    labelSize: 6
                }
            }
        )
        .highlight('cytobands', cytobands, {
            innerRadius: width/2 - 150,
            outerRadius: width/2 - 130,
            opacity: 1,
            color: function (d) {
                return gieStainColor[d.gieStain]
            },
            tooltipContent: function (d) {
                return d.name
            }
        })
        .chords(
            'l1',
            linkage,
            {
                color: function (d) {return d.color},
                logScale: false,
                zIndex: 100,
                opacity: 1,
                tooltipContent: function (d) {
                    return `${d.source.label} <-> ${d.target.label}<br> ${d.source.id}:${Math.round(d.source.position)} <-> ${d.target.id}:${Math.round(d.target.position)} <br> Jacccard: ${d.jaccard} &nbsp; Cor: ${d.cor}`
                }
            }
        )
        .scatter(
            'freqEu',
            freqEu,
            {
                color: '#969696',
                strokeColor: 'black',
                strokeWidth: 1,
                shape: 'circle',
                size: 10,
                min: 0,
                max: 1,
                innerRadius: 1.075,
                outerRadius: 1.175,
                axes: [
                    {
                        spacing: 0.5,
                        thickness: 1,
                        color: '#969696',
                        opacity: 0.3
                    }
                ],
                backgrounds: [
                    {
                        start: 0,
                        end: 1,
                        color: '#969696',
                        opacity: 0.5
                    }
                ],
                tooltipContent: function (d, i) {
                    return `${d.label} <br> ${d.block_id}:${Math.round(d.position)} <br> present percentage in <i>Eukaryota</i>: ${d.value.toFixed(2)}`
                }
            })
        .scatter(
            'freqArc',
            freqArc,
            {
                color: '#BDBDBD',
                strokeColor: 'black',
                strokeWidth: 1,
                shape: 'circle',
                size: 10,
                min: 0,
                max: 1,
                innerRadius: 1.195,
                outerRadius: 1.275,
                axes: [
                    {
                        spacing: 0.5,
                        thickness: 1,
                        color: '#BDBDBD',
                        opacity: 0.3
                    }
                ],
                backgrounds: [
                    {
                        start: 0,
                        end: 1,
                        color: '#BDBDBD',
                        opacity: 0.5
                    }
                ],
                tooltipContent: function (d, i) {
                    return `${d.label} <br> ${d.block_id}:${Math.round(d.position)} <br> present percentage in <i>Archaea</i>: ${d.value.toFixed(2)}`
                }
            })
        .scatter(
            'freqBac',
            freqBac,
            {
                color: '#D9D9D9',
                strokeColor: 'black',
                strokeWidth: 1,
                shape: 'circle',
                size: 10,
                min: 0,
                max: 1,
                innerRadius: 1.295,
                outerRadius: 1.395,
                axes: [
                    {
                        spacing: 0.5,
                        thickness: 1,
                        color: '#D9D9D9',
                        opacity: 0.3
                    }
                ],
                backgrounds: [
                    {
                        start: 0,
                        end: 1,
                        color: '#D9D9D9',
                        opacity: 0.5
                    }
                ],
                tooltipContent: function (d, i) {
                    return `${d.label} <br> ${d.block_id}:${Math.round(d.position)} <br> present percentage in <i>Bacteria</i>: ${d.value.toFixed(2)}`
                }
            })
        .render()
}

d3.queue()
    .defer(d3.json, './circosConfig/hsa_genome.json')
    .defer(d3.csv, './circosConfig/hsa_cytobands.csv')
    .defer(d3.csv, './circosConfig/linkage.csv')
    .defer(d3.csv, './circosConfig/freqEu.csv')
    .defer(d3.csv, './circosConfig/freqArc.csv')
    .defer(d3.csv, './circosConfig/freqBac.csv')
    .await(drawCircos)
