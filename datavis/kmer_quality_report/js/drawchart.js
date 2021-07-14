var visId = 'svg#vis'

var kqr = {
    'xScale': 0, // the current x scale
    'xMinValue': 0, // the smallest X value (always 0)
    'xMinScale': 0, // the value furthest to the left of the X slider
    'xMaxScale': 0, // the value furthest to the right of the X slider
    'yScale': 0, // the current y scale
    'yMinValue': 0, // the smallest allowed Y value (slightly bigger than the smallest Y)
    'yMinScale': 0, // the value furthest to the left of the Y slider
    'yMaxScale': 0, // the value furthest to the right of the Y slider
    'area': true,
    'unitFormat': 'human',
    'xDataKey': 'x',
    'xAxisLabel': 'cumulative kmer count',
    'xAxisUnits': 'human',
    'yDataKey': 'y',
    'yAxisLabel': 'count',
    'yAxisUnits': 'human',
    'targetDatapoints': 10000, // number of datapoints to shoot for when calculating bins
    'binSize': 0, // track the number of bins (autoscaled on first render)
    'maxBinSize': 10, // biggest likely bin size
    'windowSize': 100, // sliding window to use when looking for start of second peak
}

// Attempt to convert strings to a more likely type
// Only works on one-dimensional arrays of flat objects
function typeCast(ambiguous) {
    var fixt = []
    ambiguous.forEach(function(obj) {
        var o = {}
        for (var key in obj) {
            if (obj.hasOwnProperty(key)) {
                var val = obj[key]
                switch (true) {
                    case val == parseInt(val):
                        o[key] = parseInt(val)
                        break
                    case val == parseFloat(val):
                        o[key] = parseFloat(val)
                        break
                    case val === 'true':
                        o[key] = true
                        break
                    case val === 'false':
                        o[key] = false
                        break
                    default:
                        o[key] = val
                }
            }
        }
        fixt.push(o)
    })
    return fixt
}

function updateSliders() {
    $('#set_x')
        .slider('option', 'value', kqr.xScale)
        .slider('option', 'max', (kqr.xIdealScale * 5))
        .slider('option', 'min', (kqr.xMinScale))

    $('#set_y')
        .slider('option', 'value', kqr.yScale)
        .slider('option', 'max', (kqr.yIdealScale * 5))
        .slider('option', 'min', (kqr.yMinScale))

    $('#set_detail')
        .slider('option', 'value', kqr.maxBinSize - kqr.binSize + 1)
        .slider('option', 'max', kqr.maxBinSize)
        .slider('option', 'min', 1)

    $('.slider')
        .slider({
            'disabled': false
        })

    $('.scale_control')
        .prop('disabled', false)
}

function set_x_scale(scale) {
    kqr.xScale = scale

    if (scale >= 1000000) {
        kqr.xAxisUnits = 'SI'
    } else {
        kqr.xAxisUnits = 'human'
    }

    refreshChart()
}

function set_y_scale(scale) {
    kqr.yScale = scale

    if (scale >= 1000000) {
        kqr.yAxisUnits = 'SI'
    } else {
        kqr.yAxisUnits = 'human'
    }

    refreshChart()
}

// Calculate the volume under the graph for a given range.
function calcVolume(range) {
    data = kqr.original_data[0].data

    if (!range) {
        range = [0, Infinity]
    }

    if (range[1] > (data.length - 1)) {
        range[1] = data.length - 1
    }

    if (!data[range[1]].volume)
        console.warn(range[1])

    // Starting from the beginning? We already computed that.
    if (range[0] == 0) {
        return data[range[1]].volume
    }

    // We already pre-computed the running total, so the volume is
    // simply the difference between the range ends.
    return data[range[1]].volume - data[range[0]].volume
}

// Simple crossfilter.
// Returns dataset filtered by viewport constraints.
function filterData() {
    var chartdata = []

    // filter points off the right edge
    kqr.dataset.byX.filter([0, kqr.xScale + 1])

    // filter points smaller than the cutoff
    // ds.byY.filter([kqr.minimumValue, Infinity])

    // Show all (vertical) points for line graphs
    kqr.dataset.byY.filterAll()

    chartdata.push({
        key: kqr.dataset.metadata.name,
        values: kqr.dataset.byX.bottom(kqr.xScale),
        area: kqr.area
    })

    return [chartdata]
}

// For a given window_size, find the smallest volume within range of data.
function find_minimum(window_size, range, exit_early) {
    data = kqr.dataset.data

    // Can't go further than the last x value
    if (!range) {
        range = [0, (data.length - 1)]
    }

    range[1] = _.min([range[1], (data.length - 1)])

    var window_pos, sample
    var smallest = Infinity
    var best = 0
    for (window_pos = range[0]; window_pos <= range[1]; window_pos++) {
        sample = calcVolume([window_pos, window_pos + window_size])

        if (sample < smallest) {
            smallest = sample
            best = window_pos
        } else if ((exit_early) && (sample > smallest)) {
            return best
        }
    }
    return best
}

// Round up to the nearest magnitude
function roundUp(value) {
    var ne = value.toExponential(1)
        .split('e')

    return (Math.ceil(ne[0])) * Math.pow(10, parseInt(ne[1]))
}

// Attempt to find the ideal, minima, and maxima for X and Y scales
// Warning: magic and hand waving ahead!
function refreshScales() {

    // Use all datapoints for this calculation
    kqr.dataset.byX.filterAll()

    // xMinValue is always zero. Otherwise you'll get artifacts on the left side of the chart.
    kqr.xMinValue = 0

    // xMinScale should be bigger than the bin size
    kqr.xMinScale = kqr.binSize * 5

    // Find the largest possible X value
    kqr.xMax = _.max([kqr.xMax, kqr.dataset.byX.top(1)[0][kqr.xDataKey]])

    // max X scale == the largest available value of X.
    kqr.xMaxScale = _.last(kmer_data[0].data)[kqr.xDataKey]

    // To calculate the ideal X, we need to work with the graph volume.
    kqr.totalVolume = calcVolume([0, Infinity])

    // Ideal X is at 98% of the graph
    targetVolume = calcVolume([1, Infinity]) * 0.98
    for (ideal = 1; calcVolume([1, ideal]) < targetVolume; ideal++) {}

    // Set the ideal scale to something nice and round
    kqr.xIdealScale = roundUp(ideal)

    // If the current X scale isn't set, make it ideal
    if (!kqr.xScale) {
        kqr.xScale = kqr.xIdealScale
    }

    // The best minimum Y is 10 samples more than the smallest value. This handles the
    // edge case where a dataset has fewer than 10 samples.
    kqr.minSamples = 10
    kqr.minSamples = _.min([kqr.minSamples, kqr.dataset.data.length])
    kqr.yMinValue = _.max([kqr.yMinValue, kqr.dataset.byY.bottom(kqr.minSamples)[(kqr.minSamples - 1)][kqr.yDataKey]])

    // Max y scale == max y value
    kqr.yMaxScale = _.max([kqr.yMaxScale, kqr.dataset.byY.top(1)[0][kqr.yDataKey]])

    // ideal y: the largest y value rounded up
    kqr.yIdealScale = roundUp(kqr.yMaxScale)

    // Minimum useful y scale is debatable. This works pretty well.
    kqr.yMinScale = kqr.yIdealScale / 100
}

function updateCrossfilter() {

    kqr.dataset.crossfilter = crossfilter(typeCast(kqr.dataset.data))

    kqr.dataset.byX = kqr.dataset.crossfilter.dimension(function(d) {
        return d[kqr.xDataKey]
    })

    kqr.dataset.byY = kqr.dataset.crossfilter.dimension(function(d) {
        return d[kqr.yDataKey]
    })
    kqr.dataset.Total = kqr.dataset.data.length
    kqr.dataset.ready = true

    // Build a crossfilter for the original data too (for comparison to current view)
    kqr.original_data[0].crossfilter = crossfilter(typeCast(kqr.original_data[0].data))
    kqr.original_data[0].byX = kqr.original_data[0].crossfilter.dimension(function(d) {
        return d[kqr.xDataKey]
    })
    kqr.original_data[0].byY = kqr.original_data[0].crossfilter.dimension(function(d) {
        return d[kqr.yDataKey]
    })
}

// Called once at startup to create kqr.nvchart.
// If chart parameters change later, call refreshChart()
function showChart() {

    // Make sure the crossfilter is up to date
    updateCrossfilter()

    nv.addGraph(function() {
        kqr.nvchart = nv.models.lineChart()
        // kqr.nvchart = nv.models.scatterChart()
        .margin({
            left: 100,
            right: 80
        })
            .useInteractiveGuideline(true)
            .transitionDuration(50)
            .showLegend(true)
            .showYAxis(true)
            .showXAxis(true)

        kqr.nvchart.xAxis
            .axisLabel(kqr.xAxisLabel)

        kqr.nvchart.yAxis
            .axisLabel(kqr.yAxisLabel)

        d3.select(visId)
            .data(filterData())
            .call(kqr.nvchart)

        //Update the chart when window resizes.
        nv.utils.windowResize(function() {
            kqr.nvchart.update()
        })
    })

    if (kqr.dataset.ready) {
        refreshScales()

        set_x_scale(kqr.xIdealScale)
        set_y_scale(kqr.yIdealScale)

        // Give everything a quick moment to settle before refreshing.
        _.delay(function() {
            refreshChart()
        }, 100)

    }
}

function refreshChart() {

    // Only update the chart if it has been rendered.
    if (!kqr.nvchart) {
        return
    }

    kqr.datapointsShown = kqr.dataset.byX.bottom(kqr.xScale)
        .length
    kqr.datapointsAvailable = kqr.original_data[0].byX.bottom(kqr.xScale)
        .length

    kqr.nvchart.xAxisTickFormat(nv.format(kqr.xAxisUnits))
    kqr.nvchart.yAxisTickFormat(nv.format(kqr.yAxisUnits))

    kqr.nvchart
        .xDomain([kqr.xMinValue, kqr.xScale])
        .yDomain([kqr.yMinValue, kqr.yScale])

    d3.select(visId)
        .data(filterData())
        .call(kqr.nvchart)

    $('#datapoints')
        .text('Showing ' + kqr.datapointsShown + ' / ' + kqr.datapointsAvailable + ' datapoints for this range')

    updateSliders()
}

// AKA histogramify
// NOTE: This assumes that data is already sorted by xDataKey
function make_bins(bin_size, data) {
    if (bin_size < 2)
        return data

    var bin_total = 0
    var target_bin = 0
    var hist = []

    _.each(data, function(d) {
        // overshot?
        if (d[kqr.xDataKey] > target_bin) {
            // put what you've got in the current bin
            var entry = {}
            entry[kqr.xDataKey] = target_bin
            entry[kqr.yDataKey] = parseInt(bin_total / bin_size)
            hist.push(entry)
            // start the next bin
            target_bin = Math.ceil(d[kqr.xDataKey] / bin_size) * bin_size
            bin_total = d[kqr.yDataKey]
        } else {
            // accumulate until we overshoot
            bin_total = bin_total + d[kqr.yDataKey]
        }
    })

    return hist
}

function set_bin_size(bin_size) {
    // minimum bin size: 1
    var newBinSize = _.max([parseInt(bin_size), 1])

    // No change? Nothing to do.
    if (newBinSize == kqr.binSize) {
        return
    }

    // Save the current data in case we need to revert
    var olddata = $.extend(true, [], kqr.dataset.data)

    // recompute the bins
    kqr.dataset.data = make_bins(newBinSize, kqr.original_data[0].data)

    // update the crossfilter so .byX is up to date
    updateCrossfilter()
    filterData()

    if (kqr.dataset.byX.bottom(Infinity)
        .length > kqr.targetDatapoints) {
        if (!confirm("Warning: This will attempt to show " + kqr.dataset.byX.bottom(Infinity)
            .length + " datapoints at once, which may render VERY slowly. Continue?")) {
            kqr.dataset.data = olddata
            updateCrossfilter()
            filterData()
            return
        }
    }

    // Confirmed. Save the bin size and refresh the chart.
    kqr.binSize = newBinSize

    if (kqr.binSize == 1) {
        $('.fewer_bins')
            .addClass('disabled')
    } else {
        $('.fewer_bins')
            .removeClass('disabled')
    }

    refreshChart()
}

_.defer(function() {
    $('#set_x')
        .slider({
            min: kqr.xMinScale,
            max: kqr.xMaxScale,
            value: kqr.xMinScale,
            orientation: 'horizontal',
            // slide: function(event, ui) {
            //     $('#set_x_direct')
            //         .val(ui.value)
            // },
            stop: function(event, ui) {
                kqr.xScale = ui.value
                refreshChart()
            }
        })

    $('#set_y')
        .slider({
            min: kqr.yMinScale,
            max: kqr.yMaxScale,
            value: kqr.yMinScale,
            orientation: 'horizontal',
            // slide: function(event, ui) {
            //     $('#set_y_direct')
            //         .val(ui.value)
            // },
            stop: function(event, ui) {
                kqr.yScale = ui.value
                refreshChart()
            }
        })

    $('#set_detail')
        .slider({
            min: 1,
            max: kqr.maxBinSize,
            value: kqr.binSize,
            orientation: 'horizontal',
            // slide: function(event, ui) {
            //     $('#set_y_direct')
            //         .val(ui.value)
            // },
            stop: function(event, ui) {
                set_bin_size(kqr.maxBinSize - ui.value)
            }
        })

    $('.slider')
        .slider({
            'disabled': true
        })
    $('.scale_control')
        .prop('disabled', true)

    // Default visible tab: scale
    $('.tab')
        .removeClass('active')
    $('.tab-pane')
        .removeClass('active')

    $('#tab-scale')
        .addClass('active')
    $('#pane-scale')
        .addClass('active')

    // $('.hover-highlight')
    //     .tooltip()

    // kmer_data will be predefined by the server.

    // Pre-calculate the running volume of the graph.
    var runningVolume = 0

    kmer_data[0].data.forEach(function(d) {
        d.volume = d.y + runningVolume
        runningVolume = d.volume
    })

    // Make a deep copy clone for later re-binning.
    kqr.original_data = $.extend(true, [], kmer_data)
    kqr.datasets = kmer_data
    kqr.dataset = kmer_data[0]

    // Keep the number of datapoints to something sensible. This picks an initial histogram
    // bin size to keep the resulting dataset close to kqr.targetDatapoints.
    set_bin_size(parseInt(kmer_data[0].data.length / kqr.targetDatapoints))
    kqr.maxBinSize = kqr.binSize

    if (kqr.binSize == 1) {
        kqr.maxBinSize = 10
        console.log(kqr.original_data[0].data.length + ' datapoints')
    } else {
        console.log(kqr.original_data[0].data.length + ' datapoints, downsampled to ' + kqr.dataset.data.length)
    }

    // Finally, show the chart!
    showChart()

    // Reset the bin size to the best value for the current view.
    set_bin_size(parseInt(kqr.xScale / kqr.targetDatapoints))
})

//