/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.seg.segment;

import com.ciutils.cmodel.SegParams;
import com.cm.cyto.Cell;
import com.cm.cyto.CellList;
import com.cm.cyto.Chain;
import com.cm.cyto.LableCell;
import com.cm.model.AlPoint;
import com.cm.model.CIRaster;
import com.cm.model.Log;
import com.cm.util.Tools;
import java.util.ArrayList;
import java.util.logging.Logger;
import com.prim.model.Polygone;
import com.seg.model.StageCoord;
import java.io.IOException;
import java.text.DecimalFormat;

/**
 *
 * @author aharriso
 */
public class Inflect {

    private static final Logger LOG = java.util.logging.Logger.getLogger("cimg");
    private final String propFile = "c:\\oncology\\props\\inflect.prp";
    private final Log logFile = null; //new Log("c:\\oncology\\logs\\inflect.log");
    private final Log logger = null; //new Log("c:\\oncology\\logs\\inflecta.log");
    private final DecimalFormat dfa = new DecimalFormat("###.00");
    private final Tools tools = new Tools();
    private Cell bcell = null;
    private final CIRaster raster = new CIRaster();
    private StageCoord stageCoord = null;
    private SegParams segParams = null;
    private int reach = 6;
    private int segBorder = 2;
    private boolean dump = false;
    private int inflectAreaProc = 1500;
    private int minChainLength = 100;
    private int maxAngle = 170;
    private double maxRange = 1000.;
    private int minArea = 200;
    private int inflectDark = 0;
    private double inflectRatio = 0.;
    private int id = 0;
    private String dumpPath = "c:\\oncology\\test2\\";

    public Inflect(int index) {
        super();
        id = index;
        init();
    }

    private void init() {
        // System.out.println("Inflect dump is "+dump);

    }

    private CellList clean(Cell ncell) {
        int i, j;
        CellList res = new CellList(10);
        CellList ilist = new CellList(2);
        CIRaster lraster = new CIRaster(ncell.width, ncell.height);
        lraster.setBuffer(ncell.getMask());
        if (dump) {
            lraster.dumpRaster(10);
        }
        Chain chain = lraster.loadChain(1);
        lraster.cleanEdge();
        ncell.setMask(lraster.getBuffer());
        if (dump) {
            lraster.dumpRaster(11);
        }
        ArrayList<LableCell> lclist = new ArrayList<>(10);
        lraster.label(lclist);
        if (lclist.size() > 1) {
            ilist = loadCells(ncell, lclist);
            Cell cell;
            int x0, x1, y0, y1;
            if (ilist.size() > 1) {
                boolean minsized = false;
                for (j = 0; j < ilist.size(); j++) {
                    cell = ilist.get(j);
                    lraster.setBuffer(cell.getMask(), cell.width, cell.height());
                    lraster.border(2);
                    cell.setMask(lraster.getBBuffer());
                    res.add(cell);
                }
            }
        } else {
            res.add(ncell);
        }

        //System.out.println("Cleaned " + res.size());
        return res;
    }

    private void drawPositions(CIRaster raster, Chain chain) {
        Polygone poly = chain.getPoly();

        int cpos;
        int ival;
        for (int i = 0; i < poly.npoints; i++) {
            cpos = poly.ypoints[i] * raster.width() + poly.xpoints[i];
            ival = i;
            raster.getBuffer()[cpos] = ival;
        }

    }

    private void drawAnglea(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(7);
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                //System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawAngle(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.angle;
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                //System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawAngleb(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(8);
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                //System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawInflections(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(5);
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                //System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawInflectionsb(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(3) + 100;
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                // System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }
    }

    private void drawInflectionsc(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(4) + 100;
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                //System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawInversions(CIRaster raster, ArrayList<AlPoint> alist) {
        int cpos;
        int ival;
        AlPoint ivals;
        int oob = 0;
        for (int i = 0; i < alist.size(); i++) {
            ivals = alist.get(i);
            cpos = ivals.y * raster.width() + ivals.x;
            ival = ivals.get(6);
            if (ival == 0) {
                ival = 9;
            }
            if (cpos >= 0) {
                raster.getBuffer()[cpos] = ival;
            } else {
                oob++;
                // System.out.println("Found oob " + oob + " vals " + ivals.x + ", " + ivals.y);
            }
        }
        if (oob > 0) {
            //System.out.println("Found oob " + oob);
        }

    }

    private void drawCodes(CIRaster raster, Chain chain) {
        Polygone poly = chain.getPoly();
        int[] codes = chain.getICodes();
        int cpos;
        int ival;

        for (int i = 0; i < codes.length; i++) {
            cpos = poly.ypoints[i] * raster.width() + poly.xpoints[i];
            ival = codes[i];
            if (ival == 0) {
                ival = 8;
            }
            raster.getBuffer()[cpos] = ival;
        }

    }

    private void drawChain(CIRaster raster, Chain chain, byte mv) {
        Polygone poly = chain.getPoly();
        int cpos;
        int ival;

        for (int i = 0; i < chain.size(); i++) {
            cpos = poly.ypoints[i] * raster.width() + poly.xpoints[i];

            raster.getBuffer()[cpos] = mv;
        }

    }

    public float getCutMean(byte[] raw, AlPoint pa, AlPoint pb) {
        return raster.getLineMean(raw, pa.x, pa.y, pb.x, pb.y);
    }

    public StageCoord getStageCoord() {
        return stageCoord;
    }

    public CellList process(CellList nlist) {
        int i, j;
        CellList ret = new CellList();
        CellList ilist;

        for (i = 0; i < nlist.size(); i++) {
            if (nlist.get(i).getArea() > inflectAreaProc) {
                ilist = multiProcess(nlist.get(i));
                for (j = 0; j < ilist.size(); j++) {
                    ret.add(ilist.get(j));
                }
            } else {
                //nlist.get(i).setGroup(8);
                ret.add(nlist.get(i));
            }
        }
        return ret;
    }

    /**
     * This method will process divided objects until no more divisions are
     * made.
     *
     * @param newValue Cell
     * @return CellList
     */
    public CellList multiProcess(Cell newValue) {
        //save();
        CellList ret = new CellList();
        CellList ctemp;
        CellList temp = new CellList();
        CellList tempa = new CellList();
        CellList tempd = new CellList();
        CellList tempb, tempc;
        int i, j, k, n;
        int fcount = 1;
        int loops = 0;
        int cdump = 0;
        ctemp = process(newValue);
        fcount = ctemp.size();
        for (j = 0; j < ctemp.size(); j++) {
            if (ctemp.get(j).getArea() > minArea) {
                temp.add(ctemp.get(j));
            }
        }
        if (ctemp.size() > 1) {
            if (logger != null) {
                System.out.println("Found " + ctemp.size() + " collected " + temp.size());
                ctemp.saveToFileName(dumpPath + "founda" + (fcount++) + ".cmg");
                temp.saveToFileName(dumpPath + "foundb" + (fcount++) + ".cmg");
            }
        }
        ctemp.clear();

        if ((fcount == 0) || (temp.size() == 0)) {
            return ret;
        }
        tempc = new CellList(100);
        if (fcount == 1) {
            ret.add(temp.get(0)); // stop processing
            temp.clear();
        } else {
            while (temp.size() > 0) {
                if (logger != null) {
                    logger.write("Loop " + (loops) + " count " + temp.size());
                }
                if (temp.size() > 0) {
                    if (logger != null) {
                        logger.write("Collected at "+loops+"  size " + temp.size());
                    }
                    if (dump) {
                        temp.saveToFileName(dumpPath + "foundc" + (cdump++) + ".cmg");
                    }
                    loops++;
                    for (i = 0; i < temp.size(); i++) { // process all detected again
                        if (logger != null) {
                            logger.write("Cella " + i + " area " + temp.get(i).getArea() + " areaproc " + inflectAreaProc + " min " + minArea + " size " + temp.size());
                        }
                        if (temp.get(i).getArea() > inflectAreaProc) {
                            tempa = process(temp.get(i));
                            if (logger != null) {
                                logger.write("collected from " + i + " new size " + tempa.size());
                            }
                            if (dump) {
                                tempa.saveToFileName(dumpPath + "founda" + (cdump++) + ".cmg");
                            }
                            if (tempa.size() > 1) {
                                for (j = 0; j < tempa.size(); j++) {
                                    if (tempa.get(j).getArea() > minArea) {
                                        if (logger != null) {
                                            logger.write("Cellb " + j + " area " + tempa.get(j).getArea());
                                        }
                                        tempc.add(tempa.get(j)); // need to redo these
                                    } else {
                                        if (logger != null) {
                                            logger.write("Cellc " + j + " area " + tempa.get(j).getArea());
                                        }
                                    }
                                }
                            } else {
                                if (tempa.size() > 0) {
                                    if (tempa.get(0).getArea() > minArea) {
                                        ret.add(tempa.get(0)); // done with this
                                    }
                                    if (logger != null) {
                                        logger.write("Finshed " + i + " area " + tempa.get(0).getArea() + " areaproc " + inflectAreaProc + " min " + minArea);
                                    }
                                }
                            }

                        } else {
                            ret.add(temp.get(i)); // done with this
                            if (logger != null) {
                                logger.write("Finshed " + i + " area " + temp.get(i).getArea() + " areaproc " + inflectAreaProc + " min " + minArea);
                            }
                        }
                    }
                } 
                /*else {
                    if (temp.get(0).getArea() > minArea) {
                        if (logger != null) {
                            logger.write("Finished " + temp.get(0).getArea());
                        }
                        ret.add(temp.get(0)); //  single found 
                    }
                }
                */

                temp.clear();
                for (j = 0; j < tempc.size(); j++) {
                    if (tempc.get(j).getArea() > minArea) {
                        temp.add(tempc.get(j)); //  reprocess these
                        if (logger != null) {
                            logger.write("Loaded " + j + " area " + tempc.get(j).getArea() + " areaproc " + inflectAreaProc + " min " + minArea);
                        }
                    }
                }
                tempc.clear();
                tempa.clear();
                if (logger != null) {
                    logger.write("****");
                }
            }
        }
        ret.loadChains();
        return ret;
    }

    public CellList process(Cell newValue) {
        Cell cstart = newValue;
        int i;
        int j;
        boolean split = false;
        double distance;
        double minDistance;
        float meancut;
        int startArea;
        boolean inside = false;
        int chLength = minChainLength;
        CellList ilist = new CellList(2);
        ArrayList<LableCell> lclist;
        Chain chain;

        if (cstart.getArea() > inflectAreaProc) {
            bcell = new Cell(cstart);
            float cutLevel = newValue.getBackground() - inflectDark;
            raster.setSize(bcell.width, bcell.height);
            if (dump) {
                raster.setBuffer(bcell.getData(0));
                raster.dumpRaster(dumpPath + "testa.image");
            }
            raster.setBuffer(bcell.getMask(0));
            if (dump) {
                raster.dumpRaster(dumpPath + "testa.mask");
            }
            raster.setTrimBorder(0);
            raster.setSegBorder(segBorder);

            startArea = raster.getArea();
            if (dump) {
                raster.dumpRaster(1);
            }
            raster.clearEdge();
            if (dump) {
                raster.dumpRaster(2);
            }
            double circ = bcell.getCircularity();
            chain = raster.loadChain(1);
            bcell.setChain(chain);

            if ((startArea - raster.getArea()) > 0) {
                lclist = new ArrayList<>(10);
                raster.setAreaMin(20); // leave small areas
                raster.label(lclist);
                if (lclist.size() > 1) {
                    ilist = loadCells(newValue, lclist);
                    if (logger != null) {
                        logger.write("edged " + ilist.size());
                    }
                    return ilist;
                } else {
                    if (ilist.size() == 1) {
                        bcell = ilist.get(0);
                        raster.setBuffer(bcell.getIMask(0), bcell.width, bcell.height);
                        chain = raster.loadChain(1);
                    }
                }
            }
            raster.setAreaMin(minArea);
            int separation;
            int eseparation;
            if ((chain != null) && (chain.size() > 20)) {
                if (dump) {
                    chain.save(dumpPath + "prechain.chn");
                }
                for (int k = 0; k < segParams.getInflectSearch(); k++) {
                    chain.addPoint(chain.getCode(k + 1));
                }
                chain.setReach(segParams.getInflectSearch());
                chain.setMaxAngle(segParams.getInflectMaxAngle());
                chain.setInflectionLength(segParams.getInflectChainLen());

                AlPoint pa = new AlPoint();
                AlPoint pb = new AlPoint();
                AlPoint pma = new AlPoint();
                AlPoint pmb = new AlPoint();

                ArrayList<AlPoint> alist = chain.getInflections();
                pma.pos = 0;
                pmb.pos = alist.size() / 2;
                raster.reset();

                if (alist.size() > 0) {
                    if (dump) {
                        raster.dumpRaster(0);

                        drawCodes(raster, chain);
                        raster.dumpRaster(-1);
                        drawChain(raster, chain, (byte) 1); // restore border
                        drawInflections(raster, alist);
                        raster.dumpRaster(-2);
                        drawInflectionsb(raster, alist);
                        raster.dumpRaster(-3);
                        drawChain(raster, chain, (byte) 1); // restore border
                        drawInflectionsc(raster, alist);
                        raster.dumpRaster(-4);
                        drawChain(raster, chain, (byte) 1); // restore border
                        drawInversions(raster, alist);
                        raster.dumpRaster(-5);
                        drawPositions(raster, chain);
                        raster.dumpRaster(-6);
                        drawChain(raster, chain, (byte) 1); // restore border
                        drawAngle(raster, alist);
                        raster.dumpRaster(-7);
                        drawChain(raster, chain, (byte) 1); // restore border
                        drawAngleb(raster, alist);
                        raster.dumpRaster(-8);
                        drawChain(raster, chain, (byte) 1); // restore border
                    }
                    minDistance = 20000;
                    if (chLength > chain.size() / 3) {
                        chLength = chain.size() / 3;
                        //System.out.println("Reset chain size to " + chLength);
                    }
                    if (alist.size() > 1) {
                        if (dump) {
                            System.out.println("Found inflections " + alist.size() + " in size " + bcell.width + ", " + bcell.height + " Range " + maxRange + " chainlen " + chain.size() + " min " + minChainLength);
                            chain.save(dumpPath + "inchain.chn");
                        }
                        for (i = 0; i < alist.size(); i++) {
                            for (j = i + 1; j < alist.size(); j++) {
                                pa = alist.get(i);
                                pb = alist.get(j);
                                separation = pb.pos - pa.pos; // number of codes apart
                                eseparation = (chain.size() - pb.pos) + pa.pos; // wrapping
                                if (logFile != null) {
                                    //logFile.write("@ " + i + ", " + j + " sepa " + separation + " sepb " + eseparation + " min " + chLength);
                                }

                                if ((separation > chLength) && (eseparation > chLength)) {
                                    distance = Math.sqrt((pa.x - pb.x) * (pa.x - pb.x) + (pa.y - pb.y) * (pa.y - pb.y));

                                    if (dump) {
                                        System.out.println("Nodes " + pa.pos + ", " + pb.pos + " sep " + separation + " esep " + eseparation + " dist " + distance);
                                    }
                                    inside = false;
                                    if (distance < minDistance) {
                                        minDistance = distance;
                                        if (isInside(pa.x, pa.y, pb.x, pb.y)) { // check mid point is inside
                                            inside = true;
                                            meancut = 0;
                                            if (inflectDark != 0) {
                                                meancut = getCutMean(bcell.getData(0), pa, pb);
                                                if (dump) {
                                                    System.out.println("Mean cut is " + dfa.format(meancut) + " of " + dfa.format(cutLevel) + " at " + pa.pos + ", " + pb.pos + " cut " + dfa.format(distance) + " chain " + separation + " dark " + inflectDark);
                                                }
                                            }
                                            if ((inflectDark == 0) || (meancut > cutLevel)) {
                                                pma = pa;
                                                pmb = pb;
                                            }
                                        } else {
                                            if (dump) {
                                                System.out.println("Outside at " + pa.pos + " to " + pb.pos);
                                            }
                                        }
                                        if (logFile != null) {
                                            logFile.write("Inc " + pa.pos + ", " + pb.pos + " distance " + distance + " min " + minDistance + " sepa " + separation + " sepb " + eseparation + " inside " + inside);
                                        }
                                    } else {
                                        if (logFile != null) {
                                            logFile.write("rej " + pa.pos + ", " + pb.pos + " distance " + distance + " min " + minDistance + " sepa " + separation + " sepb " + eseparation);
                                        }
                                    }
                                } else {
                                    if (dump) {
                                        System.out.println(" sep " + separation + " esep " + eseparation + "chlen " + chLength);
                                    }
                                }
                            }
                        }
                        if (dump) {
                            System.out.println("Min dist " + minDistance + " Max: " + maxRange + " points " + alist.size());
                        }
                        if (minDistance < maxRange) {
                            separation = pmb.pos - pma.pos;
                            eseparation = (chain.size() - pmb.pos) + pma.pos;
                            if (logFile != null) {
                                logFile.write("last " + pma.pos + ", " + pmb.pos + " min " + minDistance + " sepa " + separation + " sepb " + eseparation + " minlength " + chLength);
                            }
                            if ((separation > chLength) && (eseparation > chLength)) {
                                split = true;
                                if (dump) {
                                    System.out.println("Cut length is " + minDistance + " minchain is " + separation);
                                    System.out.println("Drawing " + pma.x + ", " + pma.y + " to " + pmb.x + ", " + pmb.y);
                                    if (inflectDark != 0) {
                                        System.out.println("Mean cut " + getCutMean(bcell.getData(0), pma, pmb));
                                    }
                                }
                                raster.draw(pma.x, pma.y, pmb.x, pmb.y, 2); // draw blank line to separate cells
                                if (dump) {
                                    raster.dumpRaster(-9);
                                }
                                raster.draw(pma.x, pma.y, pmb.x, pmb.y, 0); // draw blank line to separate cells
                                if (dump) {
                                    raster.dumpRaster(-10);
                                }
                                if (logFile != null) {
                                    logFile.write("clip " + pma.pos + ", " + pmb.pos + " min " + minDistance + " sepa " + separation + " sepb " + eseparation + " minlength " + chLength);
                                }
                            }
                        }

                        if (split) {
                            lclist = new ArrayList<>(10);
                            raster.setAreaMin(0);
                            if (dump) {
                                raster.dumpRaster(61);
                            }
                            raster.label(lclist);
                            if (dump) {
                                System.out.println("Split found " + lclist.size() + " min area 0");
                                raster.dumpRaster(62);
                            }
                            int rcnt = lclist.size();
                            ilist = loadCells(newValue, lclist);
                            if (dump) {
                                ilist.saveToFileName(dumpPath + "splits.cmg");
                            }
                            Cell cell;
                            int x0, x1, y0, y1;
                            if (ilist.size() > 1) {
                                boolean minsized = false;
                                for (j = 0; j < 2; j++) {
                                    cell = ilist.get(j);
                                    raster.setBuffer(cell.getMask(), cell.width, cell.height());
                                    x0 = pma.x - (cell.getScreenx() - bcell.getScreenx());
                                    x1 = pmb.x - (cell.getScreenx() - bcell.getScreenx());
                                    y0 = pma.y - (cell.getScreeny() - bcell.getScreeny());
                                    y1 = pmb.y - (cell.getScreeny() - bcell.getScreeny());

                                    if (dump) {
                                        //   System.out.println("Line " + j + " @ " + x0 + ", " + y0 + ", " + x1 + ", " + y1);
                                        raster.draw(x0, y0, x1, y1, 255);
                                        raster.dumpRaster(63 + j);
                                    }
                                    raster.draw(x0, y0, x1, y1, 1); // fill dividing line
                                    raster.border(2);
                                    if (dump) {
                                        raster.dumpRaster(80 + j);
                                    }
                                    raster.border(2);
                                    cell.setMask(raster.getBBuffer());
                                    cell.setClipped(true);
                                    if (dump) {
                                        System.out.println("Area " + cell.getArea() + " min " + minArea);
                                    }
                                    if (cell.getArea() < minArea) {
                                        // minsized = true;
                                    }
                                }
                                if (minsized) {
                                    if (dump) {
                                        System.out.println("Split disabled on area min " + minArea);
                                    }
                                    ilist.clear();
                                    ilist.add(newValue);
                                }
                            }

                            //ilist.checkAreas(); // remove any blanks
                            ilist.loadChains(0);
                            if (dump) {
                                System.out.println("Split labeled " + rcnt + " found cells " + ilist.size());
                            }
                        } else {
                            ilist.add(newValue);
                        }
                    } else {
                        ilist.add(newValue);
                    }
                } else {
                    ilist.add(newValue);
                }
            } else {
                try {
                    ilist.saveCell(newValue, dumpPath + "inflect" + id + ".cmg");
                } catch (IOException io) {

                }

                //LOG.info("No inflection chain");
                ilist.add(newValue);
            }
        } else {
            if (dump) {
                System.out.println("Inflect area " + newValue.getArea() + " limit " + inflectAreaProc);
            }
            ilist.add(newValue);
        }
        if (dump) {
            //save();
        }
        ilist.loadChains();
        return ilist;
    }

    /**
     * This method will determine if the mid point is inside the object
     *
     * @param px1 int
     * @param py1 int
     * @param px2 int
     * @param py2 int
     * @return boolean
     */
    private boolean isInside(int px1, int py1, int px2, int py2) {
        boolean ret = false;
        int mid = ((py2 + py1) / 2) * raster.width() + ((px1 + px2) / 2);
        if (raster.getBuffer()[mid] != 0) {
            ret = true;
        }
        return ret;
    }

    public CellList loadCells(Cell cell, ArrayList<LableCell> alist) {
        CellList llist = new CellList(5);
        if (alist.size() > 1) {
            Cell ncell;
            LableCell lcell;
            int screenx;
            int screeny;
            int nborder;

            for (int i = 0; i < alist.size(); i++) {
                lcell = alist.get(i);
                ncell = new Cell(cell);
                ncell.setSize(lcell.width(), lcell.height());
                ncell.setNum(cell.getNum());
                //ncell.getChain().clear();
                for (int j = 0; j < cell.getImageCount(); j++) {
                    ncell.setData(j, cell.getData(j, alist.get(i).getScreenx(), lcell.getScreeny(), lcell.width(), lcell.height()));
                }
                ncell.setMask(lcell.getBMask());
                screenx = cell.getScreenx() + lcell.getScreenx();
                screeny = cell.getScreeny() + lcell.getScreeny();
                ncell.setScreen(screenx, screeny);
                if (alist.size() > 0) {
                    //ncell.setNum(12);
                }
                if ((nborder = ncell.enforceBorder(2)) > 0) {
                    //throw new NullPointerException("enforced error");
                    ncell.border(nborder);
                }
                llist.add(ncell);
            }
            if (stageCoord != null) {
                stageCoord.update(cell, llist);
            } else {
                //LOG.info("Inflect Null stageCoord");
            }
        } else {
            //cell.setNum(13);
            llist.add(cell);
        }
        llist.loadChains();
        return llist;
    }

    private void save() {
        //System.out.println("Saving " + propFile);
        tools.saveProperty(propFile, "minchain", Integer.toString(minChainLength));
        tools.saveProperty(propFile, "maxrange", Double.toString(maxRange));
        tools.saveProperty(propFile, "areaproc", Integer.toString(inflectAreaProc));
        tools.saveProperty(propFile, "minarea", Integer.toString(minArea));
        tools.saveProperty(propFile, "dark", Integer.toString(inflectDark));
        tools.saveProperty(propFile, "maxangle", Integer.toString(maxAngle));
        tools.saveProperty(propFile, "search", Integer.toString(reach));
        tools.saveProperty(propFile, "border", Integer.toString(segBorder));
    }

    public void setDump(boolean state) {
        dump = state;
    }

    public void setInfectOffset(int newValue) {
        inflectDark = newValue;
    }

    public void setMaxRange(double newValue) {
        maxRange = newValue;
    }

    public void setMinSegArea(int newValue) {
        inflectAreaProc = newValue;
    }

    public void setMinChainLen(int newValue) {
        minChainLength = newValue;
    }

    public void setSegBorder(int newValue) {
        segBorder = newValue;
    }

    public void setSegParams(SegParams newValue) {
        segParams = newValue;
        minChainLength = segParams.getInflectChainLen();
        inflectAreaProc = segParams.getInflectAreaProc();
        maxAngle = segParams.getInflectMaxAngle();
        minArea = segParams.getAreaMin();
        inflectDark = segParams.getInflectDark();
        inflectRatio = segParams.getInflectRatio();
        maxRange = segParams.getInflectRange();
        segBorder = segParams.getSegBorder();
        reach = segParams.getInflectSearch();
        // segParams.save(dumpPath+"isegs.prp");
    }

    public void setStageCood(StageCoord newValue) {
        stageCoord = newValue;
    }
}
