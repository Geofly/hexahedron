/*
 *  ::718604!
 * 
 * Copyright(C) November 20, 2014 U.S. Food and Drug Administration
 * Authors: Dr. Vahan Simonyan (1), Dr. Raja Mazumder (2), et al
 * Affiliation: Food and Drug Administration (1), George Washington University (2)
 * 
 * All rights Reserved.
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

function vjD3View ( viewer )
{
    loadCSS("d3js/css/d3js_basics.css");
    vjDataViewViewer.call(this,viewer); // inherit default behaviours of the DataViewer
    
    if(!this.downloadSvgText)this.downloadSvgText="Download graph as SVG file";
    var oThis = this;

    this.csvParserFunction=d3.csv.parse;
    this.composerFunction=function( viewer , content)
    {
        this.txt=this.getData(0).data;//content;
        this.refresh();
    };

    
    this.refresh = function ()
    {
        if (!this.div) return; //this is for the case when refresh is called from the new layout manager, sometiemes refresh is called before the div is constructed
        
        this.div.innerHTML=this.composeText(this.txt);
        //if need to use d3area.append("svg"), instead just used d3svg. not sure about this though?...
        this.d3area = d3.select("#d3js-"+this.objCls);
        this.d3svg=this.d3area.append("svg")
            .attr("style", "width='100%'; height='100%';")
            .attr("xmlns", "http://www.w3.org/2000/svg")
            .attr("xmlns:xmlns:xlink", "http://www.w3.org/1999/xlink"); // .append("svg")

        if (this.txt !== null && this.txt !== undefined) {
            if (this.newick)
                this.d3Compose(this.txt);
            else if (this.json)
                this.d3Compose(JSON.parse(this.txt));
            else if (this.csvTbl)
                this.d3Compose(this.txt);
            else if(this.d3Compose) {
                this.d3Compose(this.csvParserFunction(this.txt));
            }
            //d3.csv("?cmd=file&f=iris.csv", vjObjFunc("d3Compose",this.objCls).func );

            var that = this;
            $("#d3js-"+this.objCls+"-svg-download").on("click", function (event){
                var htmlElem = $("#"+this.id).parent().children("span").html();
                //this.href = "data:image/svg+xml," + encodeURIComponent(htmlElem);
                if (that.downloadSvg) {
                    var async = false;

                    if (!that.objectUrl) {
                        that.svgdiv = this.parentElement.children[0];
                            
                        var images = that.svgdiv.getElementsByTagNameNS(vjSVGElementSource, "image");
                        var cnt_image_fixup_needed = 0;
                        for (var i = 0; i < images.length; i++) {
                            var href = images[i].getAttributeNS("http://www.w3.org/1999/xlink", "href");
                            if (href && !href.startsWith("data:")) {
                                cnt_image_fixup_needed++;
                            }
                        }

                        if (cnt_image_fixup_needed) {
                            // need to clone the div and replace image hrefs with data uris
                            // async because chromer requires async-mode xhr for arraybuffer
                            async = true;
                            var download_text = gObject(that.container + "-downloadTextSVG");
                            if (download_text) {
                                download_text.innerHTML = that.downloadGeneratingText;
                            }

                            var cnt_image_fixup_done = 0;
                            var svgdiv = document.createElement("div");
                            svgdiv.appendChild(that.svgdiv.firstChild.cloneNode(true));
                            var images = svgdiv.getElementsByTagNameNS(vjSVGElementSource, "image");

                            for (var i = 0; i < images.length; i++) {
                                var href = images[i].getAttributeNS("http://www.w3.org/1999/xlink", "href");
                                if (href && !href.startsWith("data:")) {
                                    (function (viewer, image, href) {
                                        var xhr = new XMLHttpRequest();
                                        xhr.open("GET", href, true);
                                        xhr.responseType = "arraybuffer";
                                        xhr.onload = function() {
                                            var content_type = xhr.getResponseHeader("Content-Type") || "image/png";
                                            // ugly method to turn binary arraybuffer into base64 reasonably efficiently and avoid stack overflow
                                            var data = "data:" + content_type + ";base64," + btoa((function(u8ary) {
                                                var accum = "";
                                                for(var j = 0; j < u8ary.length; j++) {
                                                    accum += String.fromCharCode(u8ary[j]);
                                                }
                                                return accum;
                                            } (new Uint8Array(xhr.response))));
                                            image.setAttributeNS("http://www.w3.org/1999/xlink", "href", data);
                                            cnt_image_fixup_done++;

                                            if (cnt_image_fixup_done >= cnt_image_fixup_needed) {
                                                viewer.objectBlob = new Blob([svgdiv.innerHTML], {type: "image/svg+xml"});
                                                viewer.objectUrl = URL.createObjectURL(viewer.objectBlob);
                                                viewer.setdownloadSvg(svgdiv);
                                            }
                                        }
                                        xhr.send();
                                    } (that, images[i], href));
                                }
                            }
                        }
                    }

                    if (!async) {
                        that.downloadLink = this;
                        that.setdownloadSvg();
                    }
                }
            });
        }
    };

    this.composeText=function(content)
    {
        var t="<div id='d3js-div-"+this.objCls+"'>";
        t+="<span id='d3js-"+this.objCls+"'></span>";
        
        if (this.downloadSvg){
            t+="<a id='d3js-"+this.objCls+"-svg-download' href='#' download='" + sanitizeElementAttr((this.downloadFileName ? this.downloadFileName : "graphic") + ".svg") + "' style='width=100%'>" + "<img width=16 height=16 " +
                "src=./img/download.gif style='vertical-align:middle;'/>&nbsp;<small id='downloadTextSVG'>" + this.downloadSvgText + "</small></a>";                
        }
        t+= "</div>";
        
        return t;
    };

    this.setdownloadSvg = function(svgdiv) {
        if (!svgdiv) {
            svgdiv = this.svgdiv;
        }
        if (!this.objectUrl) {
            this.objectBlob = new Blob([svgdiv.innerHTML], {type: "image/svg+xml"});
            this.objectUrl = URL.createObjectURL(this.objectBlob);
        }

        var download_link = this.downloadLink;
        delete download_link.onclick;
        download_link.download = sanitizeElementAttr((this.downloadFileName ? this.downloadFileName : "graphic") + ".svg");
        download_link.href = this.objectUrl;
        if (__isIE && navigator.msSaveBlob) {
            // IE10 doesn't support download from createObjectURL(); instead, it uses a completely non-standard msSaveBlob API
            // https://msdn.microsoft.com/en-us/library/hh779016(v=vs.85).aspx
            var that = this;
            download_link.onclick = function() {
                window.navigator.msSaveBlob(that.objectBlob, download_link.download);
            }
        }
        var download_text = gObject(this.objCls + "-svg-download");
        if (download_text) {
            download_link.innerHTML = this.downloadSvgText;
        }
    };

    this.d3Compose_prv=function(flowers){}
    
    if (!d3)
        return;
        
    this.d3Tooltip = d3.select('body')
       .append('div');
       //.attr('class', 'tooltip');
    
    
     this.tooltipOn = function(d, i) {
       this.d3Tooltip.html(this.id ? d[this.id] : d ).style('visibility', 'visible');
     };
     this.tooltipMove = function(d, i) {this.d3Tooltip.style('top', (d3.event.pageY + 10) + 'px').style('left', (d3.event.pageX + 10) + 'px');}
     this.tooltipOff = function() {this.d3Tooltip.style('visibility', 'hidden');};

     

     
     return this;
}


//standard function for text wrapping
//many graphs will need this
function wrap (text, width, totalTextLength, fontOptions) {
    if (!totalTextLength)
        totalTextLength = Number.MAX_VALUE;
    if (!fontOptions){
        fontOptions = {};
    }
    if (!fontOptions.size) fontOptions.size = "12px";
    if (!fontOptions.weight) fontOptions.weight = "normal";

    text.each(function() {
        var text = d3.select(this),
            words = text.text().split(/\s+/).reverse(),
            word,
            line = [],
            lineNumber = 0,
            lineHeight = 0.4, // ems
            y = text.attr("y"),
            dy = parseFloat(text.attr("dy")),
            tspan = text.text(null).append("tspan")
                .attr("x", 0)
                .attr("y", y)
                .attr("dy", dy + "em")
                .attr("font-weight", fontOptions.weight)
                .attr("font-size", fontOptions.size);
        var tmpLength = 0;
        while (word = words.pop()) {
            tmpLength += (word + " ").length;
    
            if (tmpLength == totalTextLength)
                return;
            else if (tmpLength > totalTextLength)
            {
                var lengthNeeded = word.length - (tmpLength - totalTextLength);
                var nWord = word.substring(0, lengthNeeded);
                line.push(nWord+"...");
                tspan.text(line.join(""));
                return;
            }
            line.push(word + " ");
            tspan.text(line.join(""));
            if (tspan.node().getComputedTextLength() > width || (tspan.node().getComputedTextLength() == 0 && tmpLength*12*dy > width)) {
                line.pop();
                tspan.text(line.join(" "));
                line = [word];
                tmpLength = word.length;
                tspan = text.append("tspan")
                    .attr("x", 0)
                    .attr("y", y)
                    .attr("dy", ++lineNumber * lineHeight + dy + "em").text(word)
                    .attr("font-weight", fontOptions.weight)
                    .attr("font-size", fontOptions.size);
            }
        }
    });
};


function vjD3CategoryView ( viewer )
{
    vjD3View.call(this,viewer); // inherit default behaviours of the DataViewer


    this.categoryColors=gClrTable;//["#800", "#080;", "#008"];
    this.shapes = ['circle', 'cross', 'triangle-up', 'triangle-down', 'diamond', 'square'];
    this.colorshift=3;
    
    this.d3Compose=function(flowers) 
    {
        // read species automatically
        if(!this.id) {
            var pos=this.txt.indexOf("\n");
            if(pos!=-1)
                this.id=this.txt.substring(0,pos).split(",")[0];
        }
        if(!this.species) {
            var pos=this.txt.indexOf("\n");
            if(pos!=-1)
                this.species=this.txt.substring(0,pos).split(",")[1];
        }
        
        // Autogenerate categories
        if(!this.categories) {
            this.categories=new Array();
            for( var k=0; k<flowers.length ; ++k){
                var f= flowers[k];
                var ic;for(ic=0;ic<this.categories.length && this.categories[ic]!=f[this.species];++ic);
                if(ic<this.categories.length) continue;
                this.categories.push(f[this.species]);
            }
        }
        if(!this.traits) { 
            this.traits=new Array();
            for ( f in flowers[0] ) {
                if(f==this.species || f==this.id)continue;
                var ic;for(ic=0;ic<this.traits.length && this.traits[ic]!=f;++ic);
                this.traits .push(f);
            }
        }
    }
    
    this.clr=function(d){
        var ic;for(ic=0;ic<this.categories.length && this.categories[ic]!=d[this.species];++ic);
        var clr=this.categoryColors[(this.colorshift+ic)%this.categoryColors.length];
        return clr;
    }

    this.shape=function(d){
        return d3.svg.symbol().type('triangle-up').size(18);
    }
    
    this.d3Compose_prv=this.d3Compose;

}




//http://hsilomedus.me/wp-content/uploads/d3electricField/electricField.html
//http://is.muni.cz/www/98951/43999399/51174297/data-visualization-with-d3-js-cookbook.pdf
//http://www.nytimes.com/interactive/2012/10/15/us/politics/swing-history.html
// http://bl.ocks.org/billierinaldi/raw/3779574/


//Tooltips for d3.js SVG visualizations

/*

var tip = d3.tooltip()
    .attr('class', 'd3-tip') // define the style through css
    .offset([-10, 0])
    .html(function(d) {
        return "<strong> Blabla: </strong> <span style='color:red'>" + d + "</span><br>";
    })
    
svg.call(tip);  // svg placeholder

==> 

 YOURSVGELEMENT.on("mouseover",tip.show)
               .on("mouseout", tip.hide);
               
*/

d3.tooltip = function() {
      var direction = d3_tip_direction,
          offset    = d3_tip_offset,
          html      = d3_tip_html,
          node      = initNode(),
          svg       = null,
          point     = null,
          target    = null;
    
      function tooltip(vis) {
        svg = getSVGNode(vis);
        point = svg.createSVGPoint();
        document.body.appendChild(node);
      }
    
      // show
      tooltip.show = function() {
        var args = Array.prototype.slice.call(arguments);
        if(args[args.length - 1] instanceof SVGElement) target = args.pop();
    
        var content = html.apply(this, args),
            poffset = offset.apply(this, args),
            dir     = direction.apply(this, args),
            nodel   = d3.select(node), i = 0,
            coords;
    
        nodel.html(content)
          .style({ opacity: 1, 'pointer-events': 'all',"z-index": 9999 });
    
        while(i--) nodel.classed(directions[i], false);
        coords = direction_callbacks.get(dir).apply(this);
        nodel.classed(dir, true).style({
          top: (coords.top +  poffset[0]) + 'px',
          left: (coords.left + poffset[1]) + 'px'
        });
    
        return tooltip;
      }
    
      // hide
      tooltip.hide = function() {
        nodel = d3.select(node);
        nodel.style({ opacity: 0, 'pointer-events': 'none' });
        return tooltip;
      }
    
    
      // n - name of the attribute
      // v - value of the attribute
      tooltip.attr = function(n, v) {
        if (arguments.length < 2 && typeof n === 'string') {
          return d3.select(node).attr(n);
        } else {
          var args =  Array.prototype.slice.call(arguments);
          d3.selection.prototype.attr.apply(d3.select(node), args);
        }
    
        return tooltip;
      }
    
    
      // n - name of the property
      // v - value of the property
    
      tooltip.style = function(n, v) {
        if (arguments.length < 2 && typeof n === 'string') {
          return d3.select(node).style(n);
        } else {
          var args =  Array.prototype.slice.call(arguments);
          d3.selection.prototype.style.apply(d3.select(node), args);
        }
        return tooltip;
      }
    
      // Set or get the direction of the tooltip
      // v - One of n(north), s(south), e(east), or w(west), nw(northwest),
      //     sw(southwest), ne(northeast) or se(southeast)
      tooltip.direction = function(v) {
        if (!arguments.length) return direction;
        direction = v == null ? v : d3.functor(v);
        return tooltip;
      }
    
      // Sets or gets the offset of the tip
      // v - Array of [x, y] offset
      tooltip.offset = function(v) {
        if (!arguments.length) return offset;
        offset = v == null ? v : d3.functor(v);
        return tooltip;
      }
    
      // Public: sets or gets the html value of the tooltip
      //
      // v - String value of the tip
      //
      // Returns html value or tip
      tooltip.html = function(v) {
        if (!arguments.length) return html;
        html = v == null ? v : d3.functor(v);
    
        return tooltip;
      }
    
      function d3_tip_direction() { return 'n' ;}
      function d3_tip_offset() { return [0, 0] ;}
      function d3_tip_html() { return ' '; }
    
      var direction_callbacks = d3.map({
        n:  direction_n,
        s:  direction_s,
        e:  direction_e,
        w:  direction_w,
        nw: direction_nw,
        ne: direction_ne,
        sw: direction_sw,
        se: direction_se
      }),
    
      directions = direction_callbacks.keys();
    
      function direction_n() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.n.y - node.offsetHeight,
          left: bbox.n.x - node.offsetWidth / 2
        };
      }
    
      function direction_s() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.s.y,
          left: bbox.s.x - node.offsetWidth / 2
        };
      }
    
      function direction_e() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.e.y - node.offsetHeight / 2,
          left: bbox.e.x
        };
      }
    
      function direction_w() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.w.y - node.offsetHeight / 2,
          left: bbox.w.x - node.offsetWidth
        };
      }
    
      function direction_nw() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.nw.y - node.offsetHeight,
          left: bbox.nw.x - node.offsetWidth
        };
      }
    
      function direction_ne() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.ne.y - node.offsetHeight,
          left: bbox.ne.x
        };
      }
    
      function direction_sw() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.sw.y,
          left: bbox.sw.x - node.offsetWidth
        };
      }
    
      function direction_se() {
        var bbox = getScreenBBox();
        return {
          top:  bbox.se.y,
          left: bbox.e.x
        };
      }
    
      function initNode() {
        var node = d3.select(document.createElement('div'));
        node.style({
          position: 'absolute',
          opacity: 0,
          pointerEvents: 'none',
          boxSizing: 'border-box'
        });
    
        return node.node();
      }
    
      function getSVGNode(el) {
        el = el.node();
        if(el.tagName.toLowerCase() == 'svg')
          return el;
    
        return el.ownerSVGElement;
      }
    
      // gets the screen coordinates of a shape
      //
      // Given a shape on the screen, will return an SVGPoint for the directions
      // n(north), s(south), e(east), w(west), ne(northeast), se(southeast), nw(northwest),
      // sw(southwest).
      //
      //    +-+-+
      //    |   |
      //    +   +
      //    |   |
      //    +-+-+
      //
      // Returns an Object {n, s, e, w, nw, sw, ne, se}
      function getScreenBBox() {
        var targetel   = target || d3.event.target,
            bbox       = {},
            matrix     = targetel.getScreenCTM(),
            tbbox      = targetel.getBBox(),
            width      = tbbox.width,
            height     = tbbox.height,
            x          = tbbox.x,
            y          = tbbox.y,
            scrollTop  = document.documentElement.scrollTop || document.body.scrollTop,
            scrollLeft = document.documentElement.scrollLeft || document.body.scrollLeft;
    
    
        point.x = x + scrollLeft;
        point.y = y + scrollTop;
        bbox.nw = point.matrixTransform(matrix);
        point.x += width;
        bbox.ne = point.matrixTransform(matrix);
        point.y += height;
        bbox.se = point.matrixTransform(matrix);
        point.x -= width;
        bbox.sw = point.matrixTransform(matrix);
        point.y -= height / 2;
        bbox.w  = point.matrixTransform(matrix);
        point.x += width;
        bbox.e = point.matrixTransform(matrix);
        point.x -= width / 2;
        point.y -= height / 2;
        bbox.n = point.matrixTransform(matrix);
        point.y += height;
        bbox.s = point.matrixTransform(matrix);
    
        return bbox;
      }

  return tooltip;
};




