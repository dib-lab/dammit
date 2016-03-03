var linearTrackDefaults = {
    width: 940,
    height: 500,
    left_margin: 15,
    right_margin: 15,
    bottom_margin: 5,
    axis_height: 50,
    name: "defaultlinear",
    dragresize: true
};

function genomeTrack(layout,tracks) {

    this.tracks = tracks;
    this.layout = layout;
    this.numTracks = this.countTracks();

    if('undefined' !== typeof layout) {
	// Copy over any defaults not passed in
	// by the user
	for(var i in linearTrackDefaults) {
	    if('undefined' == typeof layout[i]) {
		this.layout[i] = linearTrackDefaults[i];
	    }
	}
    }

    if('undefined' == typeof layout.plotid) {
	this.layout.plotid = layout.container.slice(1);
    }

    this.layout.containerid = layout.container.slice(1);

    this.layout.width_without_margins =
	this.layout.width - this.layout.left_margin -
	this.layout.right_margin;

    this.layout.height_without_axis = this.layout.height -
	this.layout.axis_height;

    this.itemRects = [];

    // Start with showing the entire genome unless otherwise stated
    this.visStart = 'undefined' !== typeof layout.initStart ? layout.initStart : 0;
    this.visEnd = 'undefined' !== typeof layout.initEnd ? layout.initEnd : layout.genomesize;

    this.x = d3.scale.linear()
//	.domain([0, layout.genomesize])
	.domain([this.visStart, this.visEnd])
	.range([0,this.layout.width_without_margins]);
    this.x1 = d3.scale.linear()
	.range([0,this.layout.width_without_margins])
//       	.domain([0, layout.genomesize]);
       	.domain([this.visStart, this.visEnd]);
    this.y1 = d3.scale.linear()
	.domain([0,this.numTracks])
	.range([0,(this.layout.height_without_axis-this.layout.bottom_margin)]);
    // We need x1 and y1 in the initialization's scope too
    // to deal with passing it in to make the lollipops

    this.zoom = d3.behavior.zoom()
	.x(this.x1)
	.on("zoomstart", function () {d3.event.sourceEvent.preventDefault()} )
	.on("zoom", this.rescale.bind(this))
	.on("zoomend", this.callBrushFinished.bind(this) );

    if('undefined' == typeof layout.plotid) {
	this.layout.plotid = layout.container.slice(1);
    }

    d3.select(layout.container).select("svg").remove();

    this.chart = d3.select(layout.container)
	.append("svg")
	.attr("id", function() { return layout.container.slice(1) + "_svg"; })
	.attr("width", this.layout.width)
	.attr("height", this.layout.height)
	.attr("class", "mainTracks")
	.call(this.zoom);

    this.defs = this.chart.append("defs");
    this.clipPath = this.defs.append("clipPath")
	.attr("id", "trackClip_" + this.layout.containerid)
	.append("rect")
	.attr("width", this.layout.width_without_margins)
	//	.attr("width", this.layout.width_without_margins + this.layout.right_margin)
	.attr("height", this.layout.height)
	.attr("transform", "translate(0,0)");
    //	.attr("transform", "translate(" + this.layout.left_margin + ",0)");
    
    this.drawFeatures();

    this.main = this.chart.append("g")
       	.attr("transform", "translate(" + this.layout.left_margin + ",0)")
	.attr("width", this.layout.width_without_margins)
	.attr("height", this.layout.height)
	.attr("class", "mainTrack");

    // Resize dragger
    if(this.layout.dragresize == true) {
	var dragright = d3.behavior.drag()
	    .on("dragstart", function() {  d3.event.sourceEvent.stopPropagation(); })
	    .on("drag", this.dragresize.bind(this));

	this.dragbar_y_mid = this.layout.height/2;
	this.dragbar = this.chart.append("g")
	    .attr("transform", "translate(" + (this.layout.width-this.layout.right_margin) + "," + (this.dragbar_y_mid-10) + ")")
	    .attr("width", 25)
	    .attr("height", 20)
	    .attr("fill", "lightblue")
	    .attr("fill-opacity", .2)
	    .attr("cursor", "ew-resize")
	    .attr("id", "dragbar_" + this.layout.containerid)
	    .call(dragright);

	this.dragbar.append("rect")
	    .attr("width", 25)
	    .attr("height", 20)
	    .attr("fill-opacity", 0);

	this.dragbar.append("line")
	    .attr("x1", 6)
	    .attr("x2", 6)
	    .attr("y1", 0)
	    .attr("y2", 20)
	    .attr("class", "dragbar-line");
	this.dragbar.append("line")
	    .attr("x1", 9)
	    .attr("x2", 9)
	    .attr("y1", 0)
	    .attr("y2", 20)
	    .attr("class", "dragbar-line");
	this.dragbar.append("line")
	    .attr("x1", 12)
	    .attr("x2", 12)
	    .attr("y1", 0)
	    .attr("y2", 20)
	    .attr("class", "dragbar-line");
    }
	

    this.genomesize = layout.genomesize;

    this.tip = d3.tip()
	.attr('class', 'd3-tip')
	.offset([-10, 0])
	.html(function(d) {
		return "<strong>Name:</strong> <span style='color:red'>" + d.name + "</span>";
	    });
    
    this.chart.call(this.tip);

    this.axisContainer = this.chart.append("g")
	.attr('class', 'trackAxis')
	.attr('width', this.layout.width_without_margins)
	.attr("transform", "translate(" + (this.layout.left_margin + 5) + "," + this.layout.height_without_axis + ")");

    this.xAxis = d3.svg.axis().scale(this.x1).orient("bottom")
	.innerTickSize(-this.layout.height)
	.outerTickSize(0)
	.tickFormat(d3.format("s"));

    this.axisContainer.append("g")
	.attr("class", "xaxislinear")
	.attr('width', this.layout.width_without_margins)

	.attr("transform", "translate(0," + 10 + ")")
	.call(this.xAxis);


    for(var i=0; i < this.tracks.length; i++) {
	// We're going to see what type of tracks we have
	// and dispatch them appropriately

	 if("undefined" !== typeof this.tracks[i].skipLinear
	    &&  this.tracks[i].skipLinear == true) {
	     continue;
	 }

	 if("undefined" == typeof this.tracks[i].trackFeatures) {
	     this.tracks[i].trackFeatures = "simple";
	 } else if(this.tracks[i].trackFeatures == "complex") {
	     // We need to pre-calculate the stacking order for
	     // all arrow type features

	     // 0.77 * 0.9
	     if(this.tracks[i].trackType == "stranded") {
		 this.tracks[i].baseheight = (this.y1(1) * 0.693);
	     } else {
		 this.tracks[i].baseheight = (this.y1(1) * 0.616);
	     }
	     var inframe = [];
	     this.tracks[i].maxStackOrder = 1;
	     for(var j = 0; j < this.tracks[i].items.length; j++) {
		 // If it's the arrow type we're looking for
		 if("undefined" !== typeof this.tracks[i].items[j].feature && this.tracks[i].items[j].feature == "arrow") {
		     item = this.tracks[i].items[j];
		     // Is there anything else still in frame we need to check?

		     this.tracks[i].items[j].stackOrder = 1;
		     // Next one we encoutner will be 2
		     var curr_stackorder = 2;
		     // Go backwards through the possible inframe elements
		     // and increase their stack order if needed
		     for(var k = inframe.length - 1; k >= 0; k--) {
			 curr_k = inframe[k];
			 if(this.tracks[i].items[curr_k].end >= item.start) {
			     // If the current item in the possible stack
			     // items is overlapping...
			     if(this.tracks[i].items[curr_k].stackOrder <= curr_stackorder) {
				 // If the examined item is below or
				 // equal to the current stack order
				 this.tracks[i].items[curr_k].stackOrder = curr_stackorder;
				 curr_stackorder++;
			     }
			     // This takes care of ignoring items that
			     // are already well above the current stack order
			     
			 } else {
			     // Or if it no longer overlaps, remove it
			     inframe.splice(k, 1);
			 }
		     }

		     // Push ourselves on the inframe so the next
		     // guy can check us out
		     inframe.push(j);
		     this.tracks[i].maxStackOrder = Math.max(this.tracks[i].maxStackOrder, curr_stackorder);
		     // Save the maximum stackorder we've seen
		 }
	     }
	     this.tracks[i].increment = this.y1(0.23) / this.tracks[i].maxStackOrder;
	 }

	 if("undefined" == typeof this.tracks[i].featureThreshold) {
	     this.tracks[i].featureThreshold = this.genomesize;
	 }

	switch(this.tracks[i].trackType) {
	case "gap":
	    this.itemRects[i] = this.main.append("g")
		.attr("class", this.tracks[i].trackName)
		.attr("width", this.layout.width_without_margins)
		.attr("clip-path", "url(#trackClip_" + this.layout.containerid + ")");
	    this.displayGapTrack(this.tracks[i], i);
	    break;
	case "stranded":
	    this.itemRects[i] = this.main.append("g")
		.attr("class", this.tracks[i].trackName)
		.attr("width", this.layout.width_without_margins)
		.attr("clip-path", "url(#trackClip_" + this.layout.containerid + ")");
	    if('undefined' !== typeof this.tracks[i].linear_skipInit && this.tracks[i].linear_skipInit) {
		break;
	    }
	    this.displayStranded(this.tracks[i], i);
	    break;
	case "track":
	    this.itemRects[i] = this.main.append("g")
		.attr("class", this.tracks[i].trackName)
		.attr("width", this.layout.width_without_margins)
		.attr("clip-path", "url(#trackClip_" + this.layout.containerid + ")");
	    if('undefined' !== typeof this.tracks[i].linear_skipInit && this.tracks[i].linear_skipInit) {
		break;
	    }
	    this.displayTrack(this.tracks[i], i);
	    break;
	case "glyph":
	    if(typeof this.tracks[i].linear_invert !== 'undefined' && this.tracks[i].linear_invert == true) {
		this.tracks[i].invert = -1;
	    } else {
		this.tracks[i].invert = 1;
	    }
	    if(typeof this.tracks[i].linear_padding !== 'undefined') {
		this.tracks[i].padding = this.tracks[i].linear_padding;
	    } else {
		this.tracks[i].padding = 0;
	    }
	    this.itemRects[i] = this.main.append("g")
		.attr("class", this.tracks[i].trackName)
		.attr("width", this.layout.width_without_margins)
		.attr("clip-path", "url(#trackClip_" + this.layout.containerid + ")");
	    this.displayGlyphTrack(this.tracks[i], i);
	    break;
	case "plot":
	    this.tracks[i].g = this.itemRects[i] = this.main.append("g")
		.attr("class", this.tracks[i].trackName)
		.attr("width", this.layout.width_without_margins)
		.attr("clip-path", "url(#trackClip_" + this.layout.containerid + ")");
	    this.tracks[i].g.append("path")
		.attr("class", this.tracks[i].trackName)
		.attr("id", this.tracks[i].trackName)
		.attr("stroke-width", 1)
		.attr("fill", "none");

	    this.displayPlotTrack(this.tracks[i], i);
	    break;
	default:
	    // Do nothing for an unknown track type
	}
    }

    //    this.main.append("g").attr("transform", "matrix(0.7, 0, 0, 0.7, 0, 0)").append("use").attr("xlink:href", "#lollipop");
    //    this.main.append("g").attr("transform", "translate(100,137)").on("click", function() { console.log("click!");}).append("use").attr("xlink:href", "#lollipop_strand_pos");

}

// We can't display all track types, or some don't
// add to the stacking (ie. graph type)

genomeTrack.prototype.countTracks = function() {
    var track_count = 0;

     for(var i=0; i < this.tracks.length; i++) {

	 if("undefined" !== this.tracks[i].skipLinear
	    &&  this.tracks[i].skipLinear == true) {
	     continue;
	 }

	switch(this.tracks[i].trackType) {
	case "stranded":
	    // a linear track counts as two
	    track_count++;
	    this.tracks[i].stackNum = track_count;
	    track_count++;
	    break;
	case "track":
	    this.tracks[i].stackNum = track_count;
	    track_count++;
	    break;
	default:
	    // Do nothing for an unknown track type
	}
    }
   
    return track_count;
}

genomeTrack.prototype.displayStranded = function(track, i) {
    var visStart = this.visStart,
    visEnd = this.visEnd,
    visRange = visEnd - visStart,
    x1 = this.x1,
    y1 = this.y1;
    var cfg = this.layout;

    // Because of how the tooltip library binds to the SVG object we have to turn it
    // on or off here rather than in the .on() call, we'll redirect the calls to
    // a dummy do-nothing object if we're not showing tips in this context.
    var tip = {show: function() {}, hide: function() {} };
    if(('undefined' !== typeof track.showTooltip) && track.showTooltip) {
	tip = this.tip;
    }

    var stackNum = this.tracks[i].stackNum;
    //    console.log(visStart, visEnd);
    var visItems = track.items.filter(function(d) {
	    if(typeof d.feature !== 'undefined' && d.feature !== 'gene') {
		if(track.featureThreshold < visRange) {
		    return false;
		}
	    }
	    return d.start < visEnd && d.end > visStart;
	});

    //    console.log(track.items);

    var rects = this.itemRects[i].selectAll("g")
    .data(visItems, function(d) { return d.id; })
	.attr("transform", function(d,i) { 
	    return "translate(" + x1(d.start) + ',' +  d.yshift + ")"; });

    // Process the changed/moved rects
    rects.selectAll("rect")
    .each(function (d) { d.width = x1(d.end + 1) - x1(d.start); })
    .attr("width", function(d) {return d.width;})
    // Yes we really don't need to set the class here again
    // except to deal with the _zoom class when zooming
    // in and out
    .attr("class", function(d) {
	return track.trackName + '_' + d.suffix + ' ' + ((d.width > 5) ? (track.trackName + '_' + d.suffix + '_zoomed') : '' ) + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');});

    // Process the text for changed/moved rects
    rects.selectAll("text")
    .attr("dx", "2px")
    .attr("dy", "0.94em")
    .each(function (d) {
	    var bb = this.getBBox();
	    var slice_length = x1(d.end) - x1(d.start) - 2; // -2 to offset the dx above
	    d.visible = (slice_length > bb.width);
	})
    .attr("class", function(d) {
	return track.trackName + '_text ' + track.trackName + '_' + d.suffix + '_text ' + (d.visible ? '' : "linear_hidden " ) + ('undefined' !== typeof d.extraclass ? d.extraclass : ''); });

    this.itemRects[i].selectAll(".arrow")
    .each(function(d) {
	    d.width = x1(d.end + 1) - x1(d.start);

	    if(d.strand == -1) {
		headxtranslate = 0;
		arrowline = "m " + d.width + ",0 " + "l 0," + (track.baseheight + d.height) + " -" + d.width + ",0";
	    } else {
		headxtranslate = d.width;
		arrowline = "m 0," + track.baseheight + "l 0,-" + (track.baseheight + d.height) + " " + d.width + ",0";
	    }

	    d3.select(this).select("path")
		.attr("d", function(d) { 
			return arrowline;
		    });

	    d3.select(this).select("use")
		.attr("transform", function(d) {
			return "translate(" + headxtranslate + "," + d.headytranslate + ")";
		    });

	});

    var entering_rects = rects.enter().append("g")
    .attr("transform", function(d,i) {
	    if(d.strand == -1) {
		ystack = stackNum;
		d.suffix = 'neg';
	    } else if(d.strand == "1") {
		ystack = stackNum -1;
		d.suffix = 'pos';
	    } else {
		ystack = stackNum - 0.3;
		d.suffix = 'none';
	    }
	    var shift_gene = 0;
	    if(typeof d.feature !== 'undefined' && d.feature == "terminator") {
		ystack = ystack - 0.5;
	    } else if (track.trackFeatures == 'complex' && d.strand == "1") {
		    var shift_gene = y1(1) * .2;

	    }
	    d.yshift = y1(ystack) + 10 + shift_gene;
	    return "translate(" + x1(d.start) + ',' +  d.yshift + ")"; })
    //	    return "translate(" + x1(d.start) + ',' +  (y1(ystack) + 10 + shift_gene) + ")"; })
    .attr("id", function(d,i) { return track.trackName + '_' + d.id; })
    .attr("class", function(d) {
	    return track.trackName + '_' + d.suffix + '_group ' + (typeof d.feature === 'undefined' ? 'gene' : d.feature); })//;
	   
    //        entering_rects
    .each(function(d) {
	    d.width = x1(d.end) - x1(d.start);

	    if (typeof d.feature === 'undefined' || d.feature == "gene") {
		d3.select(this)
		.append("rect")
		.attr("class", function(d) {

			return track.trackName + '_' + d.suffix + ' ' + ((d.width > 5) ? (track.trackName + '_' + d.suffix + '_zoomed') : '') + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');})
		.attr("width", function(d) {return d.width;})
		.attr("height", function(d) {
			if(track.trackFeatures == 'complex') { 
			    var scale_factor = 0.77;
			} else {
			    var scale_factor = 1;
			}
			return (d.strand == 0 ? .4 : .9) * scale_factor * y1(1);
		    })
    .on("click", function(d,i) {
	    if (d3.event.defaultPrevented) return; // click suppressed
	    if('undefined' !== typeof track.linear_mouseclick) {
		var fn = window[track.linear_mouseclick];
		if('object' ==  typeof fn) {
		    return fn.onclick(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    } else {
		null;
	    }
	})
    .on('mouseover', function(d) { 
	    tip.show(d);
	    if('undefined' !== typeof track.linear_mouseover) {
		var fn = window[track.linear_mouseover];
		if('object' ==  typeof fn) {
		    return fn.mouseover(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	})
    .on('mouseout', function(d) { 
	    tip.hide(d);
	    if('undefined' !== typeof track.linear_mouseout) {
		var fn = window[track.linear_mouseout];
		if('object' ==  typeof fn) {
		    return fn.mouseout(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	});

	    } else if(d.feature == "terminator") {
		if(d.strand == 1) {
		    lollipop = "#lollipop_strand_pos";
		} else {
		    lollipop = "#lollipop_strand_neg";
		}
		d3.select(this).append("use").attr("xlink:href", lollipop);

	    } else if(d.feature == "arrow") {

		if(d.strand == -1) {
		    d.height = y1(0.23) - (d.stackOrder * track.increment);
		    arrowhead = "#leftarrow";
		    headxtranslate = 0;
		    d.headytranslate = track.baseheight + d.height;
		    //		    console.log(track.increment);
		    //		    console.log(d.headytranslate);
		    arrowline = "m " + d.width + ",0 " + "l 0," + (track.baseheight + d.height) + " -" + d.width + ",0";
		    //		    console.log(arrowline);
		} else {
		    d.height = d.stackOrder * track.increment;
		    arrowhead = "#rightarrow";
		    headxtranslate = d.width;
		    d.headytranslate = d.height * -1;
		    arrowline = "m 0," + track.baseheight + "l 0,-" + (track.baseheight + d.height) + " " + d.width + ",0";

		}
		arrowclass = track.trackName + '_arrow_' + d.suffix + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');
		arrowbase = d3.select(this)
		
		arrowbase.append("path")
		.attr("class", arrowclass)
		.attr("d", function(d) { 
			// 0.77 * 0.9
			
			return arrowline;
		    })
		.attr("fill-opacity", 0)
		arrowbase.append("use").attr("xlink:href", arrowhead)
		.attr("transform", function(d) {
			return "translate(" + headxtranslate + "," + d.headytranslate + ")";
		    })
		.attr("class", arrowclass);


	    } //else
	});

    if(('undefined' !== typeof track.showLabels) && typeof track.showLabels) {
	entering_rects.each(function(d) {
		if(typeof d.feature == 'undefined' || d.feature == 'gene') {
		    d3.select(this).append("text")
			.text(function(d) {return d.name;})
			.attr("dx", "2px")
			.attr("dy", "1em")
			.each(function (d) {
				var bb = this.getBBox();
				var slice_length = x1(d.end - d.start);
				d.visible = (slice_length > bb.width);
			    })
			.attr("class", function(d) {
				return track.trackName + '_text ' +  track.trackName + '_' + d.suffix + '_text ' + (d.visible ? null : "linear_hidden"  ); });
		}
	    });
    }

    rects.exit().remove();
}

genomeTrack.prototype.displayTrack = function(track, i) {
    var visStart = this.visStart,
    visEnd = this.visEnd,
    visRange = visEnd - visStart,
    x1 = this.x1,
    y1 = this.y1;
    var cfg = this.layout;

    // Because of how the tooltip library binds to the SVG object we have to turn it
    // on or off here rather than in the .on() call, we'll redirect the calls to
    // a dummy do-nothing object if we're not showing tips in this context.
    var tip = {show: function() {}, hide: function() {} };
    if(('undefined' !== typeof track.showTooltip) && track.showTooltip) {
	tip = this.tip;
    }

    var stackNum = this.tracks[i].stackNum;
    //    console.log(visStart, visEnd, visRange);
    var visItems = track.items.filter(function(d) {
	    if(typeof d.feature !== 'undefined' && d.feature !== 'gene') {
		if(track.featureThreshold < visRange) {
		    return false;
		}
	    }
	    return d.start < visEnd && d.end > visStart;}
	);

    //    console.log(track.items);

    var rects = this.itemRects[i].selectAll("g")
    .data(visItems, function(d) { return d.id; })
    .attr("transform", function(d,i) { 
	    return "translate(" + x1(d.start) + ',' + d.yshift  + ")"; 
	});


    this.itemRects[i].selectAll("rect")
    .each(function (d) { d.width = x1(d.end) - x1(d.start); })
    .attr("width", function(d) {return d.width; })
    // Yes we really don't need to set the class here again
    // except to deal with the _zoom class when zooming
    // in and out
    .attr("class", function(d) {return track.trackName + ' ' + ((d.width > 5) ? (track.trackName + '_zoomed') : '' ) + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');});

    rects.selectAll("text")
    .attr("dx", "2px")
    .attr("dy", "1em")
    .each(function (d) {
	    var bb = this.getBBox();
	    var slice_length = x1(d.end) - x1(d.start) - 2; // -2 to offset the dx above
	    d.visible = (slice_length > bb.width);
	})
    .attr("class", function(d) {return track.trackName + '_text ' + (d.visible ? null : "linear_hidden" ); });

    this.itemRects[i].selectAll(".arrow")
    .each(function(d) {
	    //	console.log(d);
	    d.width = x1(d.end + 1) - x1(d.start);

	    d3.select(this).select("path")
		.attr("d", function(d) { 
			return "m 0," + track.baseheight + "l 0,-" + (track.baseheight + d.height) + " " + d.width + ",0";
		    });

	    d3.select(this).select("use")
		.attr("transform", function(d) {
			return "translate(" + d.width + ",-" + d.headytranslate + ")";
		    });

	});

    var entering_rects = rects.enter().append("g")
    .attr("transform", function(d,i) { 
	    ystack = stackNum;
	    shift_gene = 0;
	    if(typeof d.feature !== 'undefined' && d.feature == "terminator") {
		ystack = ystack - 0.6;
	    } else if (track.trackFeatures == 'complex') {
		    var shift_gene = y1(1) * .175;

	    }

	    d.yshift = y1(ystack) + 10 + shift_gene;
	    return "translate(" + x1(d.start) + ',' + d.yshift  + ")"; 
	})
    .attr("id", function(d,i) { return track.trackName + '_' + d.id; })
    .attr("class", function(d) {return track.trackName + '_group ' + (typeof d.feature === 'undefined' ? 'gene' : d.feature); })//;

    //    entering_rects
    .each(function (d) { 
	    d.width = x1(d.end) - x1(d.start); 

	    if (typeof d.feature === 'undefined' || d.feature == "gene") {

		d3.select(this).append("rect")

		    .attr("class", function(d) {return track.trackName + ' ' + ((d.width > 5) ? (track.trackName + '_zoomed') : '' ) + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');})
		    .attr("width", function(d) {return d.width; })
		    .attr("height", function(d) {
			    if(track.trackFeatures == 'complex') {
				var scale_factor = 0.77;
			    } else {
				var scale_factor = 1;
			    }

			    return .8 * scale_factor * y1(1);
			})
		    .on("click", function(d,i) {
			    if (d3.event.defaultPrevented) return; // click suppressed
			    if('undefined' !== typeof track.linear_mouseclick) {
				var fn = window[track.linear_mouseclick];
				if('object' ==  typeof fn) {
				    return fn.onclick(track.trackName, d, cfg.plotid);
				} else if('function' == typeof fn) {
				    return fn(track.trackName, d, cfg.plotid);
				}
			    } else {
				null;
			    }
			})
		    .on('mouseover', function(d) { 
			    tip.show(d);
			    if('undefined' !== typeof track.linear_mouseover) {
				var fn = window[track.linear_mouseover];
				if('object' ==  typeof fn) {
				    return fn.mouseover(track.trackName, d, cfg.plotid);
				} else if('function' == typeof fn) {
				    return fn(track.trackName, d, cfg.plotid);
				}
			    }	
			})
		    .on('mouseout', function(d) { 
			    tip.hide(d);
			    if('undefined' !== typeof track.linear_mouseout) {
				var fn = window[track.linear_mouseout];
				if('object' ==  typeof fn) {
				    return fn.mouseout(track.trackName, d, cfg.plotid);
				} else if('function' == typeof fn) {
				    return fn(track.trackName, d, cfg.plotid);
				}
			    }	
			})
		    
		    if('undefined' !== typeof track.showLabels) {
			//		entering_rects
			d3.select(this).append("text")
			    .text(function(d) {return d.name;})
			    .attr("dx", "2px")
			    .attr("dy", "1em")
			    .each(function (d) {
				    var bb = this.getBBox();
				    var slice_length = x1(d.end) - x1(d.start) -2 ; // -2 to offset the dx above
				    d.visible = (slice_length > bb.width);
				})
			    .attr("class", function(d) {return track.trackName + '_text ' + (d.visible ? null : "linear_hidden" ); });
		    }
	    } else if(d.feature == "terminator") {
		d3.select(this).append("use").attr("xlink:href", "#lollipop");
	    } else if(d.feature == "arrow") {

		d.height = d.stackOrder * track.increment;
		d.headytranslate = d.height;

		arrowclass = track.trackName + '_arrow' + ' ' + ('undefined' !== typeof d.extraclass ? d.extraclass : '');
		arrowbase = d3.select(this)
		
		arrowbase.append("path")
		.attr("class", arrowclass)
		.attr("d", function(d) { 
			// 0.77 * 0.9
			return "m 0," + track.baseheight + "l 0,-" + (track.baseheight + d.height) + " " + d.width + ",0";
			
		    })
		.attr("fill-opacity", 0)
		arrowbase.append("use").attr("xlink:href", "#rightarrow")
		.attr("transform", function(d) {
			return "translate(" + d.width + ",-" + d.headytranslate + ")";
		    })
		.attr("class", arrowclass);


	    } //else
	});


    rects.exit().remove();
}

genomeTrack.prototype.displayPlotTrack = function(track, i) {
    var visStart = this.visStart,
    visEnd = this.visEnd,
    x1 = this.x1,
    y1 = this.y1;

    if((typeof track.visible == 'undefined') || (track.visible == false)) {
    	return;
    }

    if((typeof track.linear_plot_width == 'undefined') || (typeof track.linear_plot_height == 'undefined')) {
//    if((typeof track.linear_plot_width == 'undefined')) {
	return;
    }

    var startItem = parseInt(visStart / track.bp_per_element);
    var endItem = Math.min(parseInt(visEnd / track.bp_per_element), track.items.length);
    var offset = ((startItem+1) * track.bp_per_element) - visStart;

    var items = track.items.filter(function(d, i) { return i >= startItem && i <= endItem } );

    track.plotScale = d3.scale.linear()
	.domain([track.plot_min, track.plot_max])
	.range([track.linear_plot_height+(track.linear_plot_width/2), track.linear_plot_height-(track.linear_plot_width/2)]);

    var lineFunction = d3.svg.line()
	.x(function(d, i) { return x1((i*track.bp_per_element)); } )
	.y(function(d, i) { return track.plotScale(d); } )
	.interpolate("linear");

     var plot = this.itemRects[i].selectAll("path")
	.attr("d", lineFunction(track.items))

//    plot.exit().remove();

}

genomeTrack.prototype.displayGapTrack = function(track, i) {
    var visStart = this.visStart,
    visEnd = this.visEnd,
    x1 = this.x1,
    y1 = this.y1;
    var cfg = this.layout;
    var self = this;

    // Because of how the tooltip library binds to the SVG object we have to turn it
    // on or off here rather than in the .on() call, we'll redirect the calls to
    // a dummy do-nothing object if we're not showing tips in this context.
    var tip = {show: function() {}, hide: function() {} };
    if(('undefined' !== typeof track.showTooltip) && track.showTooltip) {
	tip = this.tip;
    }

    if((typeof track.visible !== 'undefined') && (track.visible == false)) {
    	return;
    }

//    var gap_range = d3.range(0, this.layout.height, 5);
    var gap_range = d3.range(0, y1(this.numTracks)+(this.numTracks*3), 5);

    var items = track.items.filter(function(d) {return (d.start <= visEnd && d.start >= visStart) || (d.end <= visEnd && d.end >= visStart);});

    var gaps = this.itemRects[i].selectAll("path")
    .data(items, function(d) { return d.id; })
    .attr("transform", function(d,i) { 
	    return "translate(" + x1(d.start) + ', 0)'; 
	})
    .each(function (d) { d.width = x1(d.end) - x1(d.start); })
    .attr('d', function(d) { return self.jaggedPathGenerator(d.width, gap_range); } );

    var entering_gaps = gaps.enter().append("path")
    .each(function (d) { d.width = x1(d.end) - x1(d.start); })
    .attr('d', function(d) { return self.jaggedPathGenerator(d.width, gap_range); } )
    .attr("transform", function(d,i) { 
	    return "translate(" + x1(d.start) + ', 0)'; 
	})
    .attr("id", function(d,i) { return track.trackName + '_' + d.id; })
    .attr("class", function(d) {return track.trackName + ' linearplot ' + (typeof d.feature === 'undefined' ? 'gene' : d.feature); })//;
    .on("click", function(d,i) {
	    if (d3.event.defaultPrevented) return; // click suppressed
	    if('undefined' !== typeof track.linear_mouseclick) {
		var fn = window[track.linear_mouseclick];
		if('object' ==  typeof fn) {
		    return fn.onclick(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    } else {
		null;
	    }
	})
    .on('mouseover', function(d) { 
	    tip.show(d);
	    if('undefined' !== typeof track.linear_mouseover) {
		var fn = window[track.linear_mouseover];
		if('object' ==  typeof fn) {
		    return fn.mouseover(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	})
    .on('mouseout', function(d) { 
	    tip.hide(d);
	    if('undefined' !== typeof track.linear_mouseout) {
		var fn = window[track.linear_mouseout];
		if('object' ==  typeof fn) {
		    return fn.mouseout(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	});

    gaps.exit().remove();

}

genomeTrack.prototype.jaggedPathGenerator = function(width, data) {

    var down = [];
    //    var up = [];
    for(var i = 0; i < data.length; i++) {
	var offset = ((i % 2 === 0) ? 3 : -3);
	down.push({ x: offset, y: data[i] });
	down.unshift({ x: offset+width, y: data[i] });
    }
    down.push(down[0]);

    var generator = d3.svg.line()
    .x(function(d,i) { return d.x; })
    .y(function(d,i) { return d.y; })
    .interpolate("linear");

    return generator(down);
    //    return generator(down.concat(up));
}

genomeTrack.prototype.displayGlyphTrack = function(track, i) {
    var visStart = this.visStart,
    visEnd = this.visEnd,
    x1 = this.x1,
    y1 = this.y1;
    var cfg = this.layout;

    if((typeof track.visible !== 'undefined') && (track.visible == false)) {
    	return;
    }

    // Because of how the tooltip library binds to the SVG object we have to turn it
    // on or off here rather than in the .on() call, we'll redirect the calls to
    // a dummy do-nothing object if we're not showing tips in this context.
    var tip = {show: function() {}, hide: function() {} };
    if(('undefined' !== typeof track.showTooltip) && track.showTooltip) {
	tip = this.tip;
    }

    var items = track.items.filter(function(d) {return d.bp <= visEnd && d.bp >= visStart;});

    // When we move we need to recalculate the stacking order
    var stackCount = 0;
    for(var j = 0; j < items.length; j++) {
	if(items[j].bp < visStart || items[j].bp > visEnd) {
	    continue;
	}
	if(j < 1) {
	    items[j].stackCount = 0;
	    continue;
	}

	var dist = x1(items[j].bp) - x1(items[j-1].bp);

	if(dist < track.linear_pixel_spacing) {
	    items[j].stackCount = items[j-1].stackCount + 1;
	    continue;
	}

	items[j].stackCount = 0;
    }

    // Because SVG coordinates are from the top-left, the "height" is pixels DOWN from
    // the top of the image to start stacking the glyphs

    var glyphs = this.itemRects[i].selectAll("path")
    .data(items, function(d) { return d.id; })
    .attr("transform", function(d,i) { return "translate(" + (x1(d.bp) + track.padding) + ',' + (track.linear_height - (track.linear_glyph_buffer * d.stackCount * track.invert))  + ")"; });
   
    var entering_glyphs = glyphs.enter()
    .append('path')
    .attr('id', function(d,i) { return track.trackName + "_glyph" + d.id; })
    .attr('class', function(d) {return track.trackName + '_' + d.type + " linear_" + track.trackName + '_' + d.type; })
    .attr("d", d3.svg.symbol().type(track.glyphType).size(track.linear_glyphSize))
    .attr("transform", function(d,i) {  return "translate(" + (x1(d.bp) + track.padding) + ',' + (track.linear_height - (track.linear_glyph_buffer * d.stackCount * track.invert))  + ")"; })
    .on("click", function(d,i) {
	    if('undefined' !== typeof track.linear_mouseclick) {
		var fn = window[track.linear_mouseclick];
		if('object' ==  typeof fn) {
		    return fn.onclick(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    } else {
		null;
	    }
	})
    .on('mouseover', function(d) { 
	tip.show(d);
	    if('undefined' !== typeof track.linear_mouseover) {
		var fn = window[track.linear_mouseover];
		if('object' ==  typeof fn) {
		    return fn.mouseover(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	})
    .on('mouseout', function(d) { 
	tip.hide(d);
	    if('undefined' !== typeof track.linear_mouseout) {
		var fn = window[track.linear_mouseout];
		if('object' ==  typeof fn) {
		    return fn.mouseout(track.trackName, d, cfg.plotid);
		} else if('function' == typeof fn) {
		    return fn(track.trackName, d, cfg.plotid);
		}
	    }	
	});

    glyphs.exit()
    .remove();
    
}

genomeTrack.prototype.displayAxis = function() {
    this.axisContainer.select(".xaxislinear").call(this.xAxis);
}

    genomeTrack.prototype.update = function(startbp, endbp, params) {
    //    console.log(startbp, endbp);

    this.visStart = startbp;
    this.visEnd = endbp;

    this.zoom.x(this.x1.domain([startbp,endbp]));

    this.redraw();
}

genomeTrack.prototype.update_finished = function(startbp, endbp, params) {
    //    console.log("Thank you, got: " + startbp, endbp);

}

genomeTrack.prototype.resize = function(newWidth) {
    this.layout.width = newWidth;

    this.dragbar
    .attr("transform", "translate(" + (newWidth -
				       this.layout.right_margin) + "," + (this.dragbar_y_mid-15) + ")")

    this.layout.width_without_margins =
	this.layout.width - this.layout.left_margin -
	this.layout.right_margin;

    this.x
	.range([0,this.layout.width_without_margins]);
    this.x1
	.range([0,this.layout.width_without_margins]);

    this.chart
	.attr("width", this.layout.width)

    this.clipPath
	.attr("width", this.layout.width_without_margins)
    
    this.main
	.attr("width", this.layout.width_without_margins)

    this.redraw();

}

genomeTrack.prototype.dragresize = function(d) {
    var newWidth = d3.event.x;

    //    console.log(this.layout.containerid);
    this.resize(newWidth);
    //    d3.event.preventDefault();

}

genomeTrack.prototype.redraw = function() {

    for(var i = 0; i < this.tracks.length; i++) {

	 if("undefined" !== this.tracks[i].skipLinear
	    &&  this.tracks[i].skipLinear == true) {
	     continue;
	 }

	switch(this.tracks[i].trackType) {
	case 'gap':
	    this.displayGapTrack(this.tracks[i], i);
            break;
	case "stranded":
	    this.displayStranded(this.tracks[i], i);
	    break;
	case "track":
	    this.displayTrack(this.tracks[i], i);
	    break;
	case "glyph":
	    this.displayGlyphTrack(this.tracks[i], i);
	    break;
	case "plot":
	    this.displayPlotTrack(this.tracks[i], i);
	    break;
	default:
	    // Do nothing for an unknown track type
	}
    }
    
    this.axisContainer.select(".xaxislinear").call(this.xAxis);

}

genomeTrack.prototype.rescale = function() {
    var cfg = this.layout;

    var reset_s = 0;
    if ((this.x1.domain()[1] - this.x1.domain()[0]) >= (this.genomesize - 0)) {
	this.zoom.x(this.x1.domain([0, this.genomesize]));
	reset_s = 1;
    }

    if (reset_s == 1) { // Both axes are full resolution. Reset.
	this.zoom.scale(1);
	this.zoom.translate([0,0]);
    }
    else {
	if (this.x1.domain()[0] < 0) {
	    this.x1.domain([0, this.x1.domain()[1] - this.x1.domain()[0] + 0]);
	}
	if (this.x1.domain()[1] > this.genomesize) {
	    var xdom0 = this.x1.domain()[0] - this.x1.domain()[1] + this.genomesize;
	    this.x1.domain([xdom0, this.genomesize]);
	}
    }

    var cur_domain = this.x1.domain();
    this.visStart = cur_domain[0];
    this.visEnd = cur_domain[1];

    if('undefined' !== typeof this.callbackObj) {
	if( Object.prototype.toString.call( this.callbackObj ) === '[object Array]' ) { 
	    for(var obj in this.callbackObj) {
		if(this.callbackObj.hasOwnProperty(obj)) {
		    this.callbackObj[obj].update(this.x1.domain()[0], this.x1.domain()[1], { plotid: cfg.plotid } );
		}
	    }
	} else {
	    this.callbackObj.update(this.x1.domain()[0], this.x1.domain()[1], { plotid: cfg.plotid } );
	}
    }

    this.redraw();

}

genomeTrack.prototype.addBrushCallback = function(obj) {

    // We allow multiple brushes to be associated with a linear plot, if we have
    // a brush already, add this new one on.  Otherwise just remember it.

    if('undefined' !== typeof this.callbackObj) {

	if( Object.prototype.toString.call( obj ) === '[object Array]' ) { 
	    this.callbackObj.push(obj);
	} else {
	    var tmpobj = this.callbackObj;
	    this.callbackObj = [tmpobj, obj];
	}
    } else {
	this.callbackObj = obj;
    }

    // And make sure our new brush is updated to reflect
    // the current visible area
    obj.update(this.visStart, this.visEnd);
}

genomeTrack.prototype.callBrushFinished = function() {
    var cfg = this.layout;

    if('undefined' !== typeof this.callbackObj) {
	if( Object.prototype.toString.call( this.callbackObj ) === '[object Array]' ) { 
	    for(var obj in this.callbackObj) {
		if(this.callbackObj.hasOwnProperty(obj)) {
		    this.callbackObj[obj].update_finished(this.x1.domain()[0], this.x1.domain()[1], { plotid: cfg.plotid } );
		}
	    }
	} else {
	    this.callbackObj.update_finished(this.x1.domain()[0], this.x1.domain()[1], { plotid: cfg.plotid } );
	}
    }

}

////////////////////////////////////////////////
//
// Export functionality
//
////////////////////////////////////////////////

// Allowed export formats are 'png' and 'svg'
// The extension will be added, just summply the base
// filename.

// Saving to raster format is dependent on FileSaver.js
// and canvg.js (which include rgbcolor.js & StackBlur.js),
// they must be loaded before circularplot.js

genomeTrack.prototype.savePlot = function(scaling, filename, stylesheetfile, format) {
        // First lets get the stylesheet
    var sheetlength = stylesheetfile.length;
    var style = document.createElementNS("http://www.w3.org/1999/xhtml", "style");
	style.textContent += "<![CDATA[\n";
    for (var i=0;i<document.styleSheets.length; i++) {
	str = document.styleSheets[i].href;
	if(null == str) continue;

	if (str.substr(str.length-sheetlength)==stylesheetfile){
            var rules;
            if(document.styleSheets[i].cssRules) {
                rules = document.styleSheets[i].cssRules;
            } else if (document.styleSheets[i].rules) {
                rules = document.styleSheets[i].rules;
            }
            if(rules) {
                for (var j=0; j<rules.length;j++){
                    style.textContent += (rules[j].cssText + "\n");
                }
            }
            break;
    	}
    }
    style.textContent += "]]>";

    // Now we clone the SVG element, resize and scale it up
    var container = this.layout.container.slice(1);
    var containertag = document.getElementById(container);
    var clonedSVG = containertag.cloneNode(true);
    var svg = clonedSVG.getElementsByTagName("svg")[0];

    // Remove any hidden elements such as text that's not being shown
    var tags = svg.getElementsByClassName("linear_hidden")
    for(var i=tags.length-1; i>=0; i--) {
//	if(tags[i].getAttributeNS(null, "name") === name) {
	    tags[i].parentNode.removeChild(tags[i]);
//        }
    }

    // We need to resize the svg with the new canvas size
    svg.removeAttribute('width');
    svg.removeAttribute('height');
    svg.setAttribute('width', this.layout.width*scaling);
    svg.setAttribute('height', this.layout.height*scaling);

    // Update first g tag with the scaling
    g = svg.getElementsByTagName("g")[0];
    transform = g.getAttribute("transform");
//    g.setAttribute("transform", transform + " scale(" + scaling + ")");

    // Append the stylehsheet to the cloned svg element
    // so when we export it the style are inline and 
    // get rendered
    svg.getElementsByTagName("defs")[0].appendChild(style);

    // Fetch the actual SVG tag and convert it to a canvas
    // element
    var content = clonedSVG.innerHTML.trim();

    if(format == 'svg') {
	var a = document.createElement('a');
	a.href = "data:application/octet-stream;base64;attachment," + btoa(content);
	a.download = filename + ".svg";
	a.click();

    } else if(format == 'png') {
	var canvas = document.createElement('canvas');
	canvg(canvas, content);

	// Convert the canvas to a data url (this could
	// be displayed inline by inserting it in to an 
	// <img> tag in the src attribute, ie
	// <img src="'+imgData+'">
	var theImage = canvas.toDataURL('image/png');

	// Convert to a blob
	var blob = this.dataURLtoBlob(theImage);

	// Prompt to save
	saveAs(blob, filename);
    }
}

genomeTrack.prototype.dataURLtoBlob = function(dataURL) {
  // Decode the dataURL    
  var binary = atob(dataURL.split(',')[1]);
  // Create 8-bit unsigned array
  var array = [];
  for(var i = 0; i < binary.length; i++) {
      array.push(binary.charCodeAt(i));
  }
  // Return our Blob object
  return new Blob([new Uint8Array(array)], {type: 'image/png'});
}

genomeTrack.prototype.drawFeatures = function() {
    var x1 = this.x1;
    var y1 = this.y1;

    // Quick debugging variable to show click rects
    var opacity = 0;

    // Lollipop for terminator glyph (stranded, positive)
    var lollipop_strand_pos = this.defs.append("g").attr("id", "lollipop_strand_pos");
    lollipop_strand_pos.append("path")
	.attr("class", "lollipophead")
	.attr("d", function() {
		arc = "m 0,";
		arc += y1(1) * 0.69;
		arc += " a ";
		arc += y1(1) * 0.13;
		arc += ",";
		arc += y1(1) * 0.13;
		arc += " 0 1 1 ";
		arc += y1(1) * 0.09;
		arc += ",0";
		//		console.log(arc);
		return arc;
		//		"m 0,60 a 12,12 0 1 1 8,0"
	    });
    lollipop_strand_pos.append("path")
	.attr("class", "lollipopstemend")
	.attr("d", function() {
		line = "m " + (y1(1) * 0.09) + ",";
		line += y1(1) * 0.69;
		line += " l 0,";
		line += y1(1) * 0.70;
		//		console.log("line " + line);
		return line;
		//		"m 8,60 l 0,65"
		    });
    lollipop_strand_pos.append("path")
	.attr("class", "lollipopstemstart")
	.attr("d", function() {
		line = "m 0,";
		line += y1(1) * 0.69;
		line += " l 0,";
		line += y1(1) * 0.70;
		//		console.log(line);
		return line;
		//		"m 0,60 l 0,65"
	    });
    lollipop_strand_pos.append("rect")
	.attr("transform", function () {
		return "translate(-" + (y1(1) * 0.1) + "," + (y1(1) * 0.43) + ")";
	    })
	.attr("width", function() { 
		return y1(1) * 0.3;
	    })
	.attr("height", function() { 
		return y1(1) * 0.26;
	    })
	.attr("fill-opacity", opacity);
    lollipop_strand_pos.append("rect")
        .attr("transform", function() {
		return "translate(0, " + (y1(1) * 0.68) + ")";
		//		"translate(0,60)"
		    })
	.attr("width", 8)
	.attr("height", function() {
		return  y1(1) * 0.72;
	    })
	.attr("fill-opacity", opacity);

    // Lollipop for terminator glyph (stranded, negative)
    var lollipop_strand_neg = this.defs.append("g").attr("id", "lollipop_strand_neg");
    lollipop_strand_neg.append("path")
	.attr("class", "lollipophead")
	.attr("d", function() {
		arc = "m 0,";
		arc += y1(1) * 1.20;
		arc += " a ";
		arc += y1(1) * 0.13;
		arc += ",";
		arc += y1(1) * 0.13;
		arc += " 0 1 0 ";
		arc += y1(1) * 0.09;
		arc += ",0";
		//		console.log(arc);
		return arc;
		//		"m 0,60 a 12,12 0 1 1 8,0"
	    });
    lollipop_strand_neg.append("path")
	.attr("class", "lollipopstemstart")
	.attr("d", function() {
		line = "m " + (y1(1) * 0.09) + ",";
		line += y1(1) * 1.20;
		line += " l 0,";
		line += y1(1) * 0.71 * -1;
		//		console.log(line);
		return line;
		//		"m 8,60 l 0,65"
		    });
    lollipop_strand_neg.append("path")
	.attr("class", "lollipopstemend")
	.attr("d", function() {
		line = "m 0,";
		line += y1(1) * 1.20;
		line += " l 0,";
		line += y1(1) * 0.71 * -1;
		//		console.log(line);
		return line;
		//		"m 0,60 l 0,65"
	    });
    lollipop_strand_neg.append("rect")
	.attr("transform", function () {
		return "translate(-" + (y1(1) * 0.1) + "," + (y1(1) * 1.20) + ")";
	    })
	.attr("width", function() { 
		return y1(1) * 0.3;
	    })
	.attr("height", function() { 
		return y1(1) * 0.26;
	    })
	.attr("fill-opacity", opacity);
    lollipop_strand_neg.append("rect")
	.attr("transform", function() {
		return "translate(0, " + (y1(1) * 0.50) + ")";
		//		"translate(0,60)"
		    })
	.attr("width", 8)
	.attr("height", function() {
		return  y1(1) * 0.71;
	    })
	.attr("fill-opacity", opacity);

    // Lollipop for terminator glyph (unstranded)
    var lollipop = this.defs.append("g").attr("id", "lollipop");
    lollipop.append("path")
        .attr("class", "lollipophead")
	.attr("d", function() {
		arc = "m 0,";
		arc += y1(1) * 0.77;
		arc += " a ";
		arc += y1(1) * 0.13;
		arc += ",";
		arc += y1(1) * 0.13;
		arc += " 0 1 1 ";
		arc += y1(1) * 0.09;
		arc += ",0";
		//		console.log(arc);
		return arc;
		//		"m 0,60 a 12,12 0 1 1 8,0"
	    });
    lollipop.append("path")
	.attr("class", "lollipopstemend")
	.attr("d", function() {
		line = "m " + (y1(1) * 0.09) + ",";
		line += y1(1) * 0.77;
		line += " l 0,";
		line += y1(1) * 0.62;
		//		console.log(line);
		return line;
		//		"m 8,60 l 0,65"
		    });
    lollipop.append("path")
	.attr("class", "lollipopstemstart")
	.attr("d", function() {
		line = "m 0,";
		line += y1(1) * 0.77;
		line += " l 0,";
		line += y1(1) * 0.62;
		//		console.log(line);
		return line;
		//		"m 0,60 l 0,65"
	    });
    lollipop.append("rect")
	.attr("transform", function () {
		return "translate(-" + (y1(1) * 0.1) + "," + (y1(1) * 0.51) + ")";
	    })
	.attr("width", function() { 
		return y1(1) * 0.3;
	    })
	.attr("height", function() { 
		return y1(1) * 0.26;
	    })
	.attr("fill-opacity", opacity);
    lollipop.append("rect")
        .attr("transform", function() {
		return "translate(0, " + (y1(1) * 0.76) + ")";
		//		"translate(0,60)"
		    })
	.attr("width", 8)
	.attr("height", function() {
		return  y1(1) * 0.64;
	    })
	.attr("fill-opacity", opacity);

    // Right arrow
    var rightarrow = this.defs.append("g").attr("id", "rightarrow");
    rightarrow.append("path")
        .attr("class", "rightarrow")
        .attr("d", function() {
		return "m 0,0 l -4,4 m 0,-8 l 4,4";
	    });

    // Left arrow
    var leftarrow = this.defs.append("g").attr("id", "leftarrow");
    leftarrow.append("path")
        .attr("class", "rightarrow")
        .attr("d", function() {
		return "m 0,0 l 4,4 m 0,-8 l -4,4";
		//		return "m 0," + y1(1) + " l 4,4 m 0,-8 l -4,4";
	    });

}
