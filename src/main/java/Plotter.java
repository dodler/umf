
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.*;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;
import java.util.Collection;
import java.util.LinkedList;

/**
 * Created by artem on 18.12.15.
 */
public class Plotter extends ApplicationFrame
{
    private int width = 600, height = 800;

    private Collection<TextTitle> legend;

    public Plotter(String title)
    {
        super(title);
        legend = new LinkedList<>();
    }

    public Plotter(String title, int width, int height)
    {
        this(title);
        this.width = width;
        this.height = height;
    }

    final XYDataset dataset = new XYSeriesCollection();

    public void addRealVector(RealVector arguments, RealVector function, String label)
    {
        final XYSeries data = new XYSeries(label);
        double[] args = arguments.toArray(), foo = function.toArray();

        for (int i = 0; i < args.length; i++)
        {
            data.add(args[i], foo[i]);
        }
        ((XYSeriesCollection) dataset).addSeries(data);
        final TextTitle title = new TextTitle(label);
        title.setPosition(RectangleEdge.BOTTOM);
// legend.add(title);
    }

    public void plot()
    {
        final JFreeChart chart = createChart(dataset);
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(width, height));
        setContentPane(chartPanel);
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }

    /**
     * Creates a chart.
     *
     * @param dataset the data for the chart.
     * @return a chart.
     */
    private JFreeChart createChart(final XYDataset dataset)
    {

// create the chart...
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "График функции u(r,t)", // chart title
                "X", // x axis label
                "Y", // y axis label
                dataset, // data
                PlotOrientation.VERTICAL,
                true, // include legend
                true, // tooltips
                false // urls
        );

// NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
        chart.setBackgroundPaint(Color.white);

// final StandardLegend legend = (StandardLegend) chart.getLegend();
// legend.setDisplaySeriesShapes(true);

// get a reference to the plotAll for further customisation...
        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.lightGray);
// plotAll.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.white);

        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesLinesVisible(1, false);
        renderer.setSeriesShapesVisible(0, false);
        plot.setRenderer(renderer);

// change the auto tick unit selection to integer units only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
// OPTIONAL CUSTOMISATION COMPLETED.

        for (TextTitle tt : legend)
        {
            chart.addSubtitle(tt);
        }

        return chart;

    }
}