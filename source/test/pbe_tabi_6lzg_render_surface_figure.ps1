param(
    [string]$ComparisonDir = 'test\tmp\pbe_tabi_6lzg_electrostatic_surface_contract\comparison',
    [string]$PymolExe = 'pymol.exe',
    [string]$OutputPath = '',
    [int]$Width = 2800,
    [int]$Height = 1600
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$ComparisonDir = [System.IO.Path]::GetFullPath($ComparisonDir)
$loader = Join-Path $ComparisonDir 'view_6lzg_titan_vs_apbs.py'
$summaryPath = Join-Path $ComparisonDir '6lzg_surface_comparison_summary.csv'
if (-not (Test-Path -LiteralPath $loader)) {
    throw "Missing PyMOL loader: $loader"
}
if (-not (Test-Path -LiteralPath $summaryPath)) {
    throw "Missing surface comparison summary: $summaryPath"
}
if (-not $OutputPath) {
    $OutputPath = Join-Path $ComparisonDir '6lzg_titan_vs_apbs_surface_figure.png'
} else {
    $OutputPath = [System.IO.Path]::GetFullPath($OutputPath)
}

$resolvedPymol = Get-Command $PymolExe -ErrorAction SilentlyContinue
if (-not $resolvedPymol) {
    throw "Unable to resolve PyMOL executable: $PymolExe"
}

$rawName = '6lzg_surface_preview_raw.png'
$rawPath = Join-Path $ComparisonDir $rawName
Push-Location $ComparisonDir
try {
    & $resolvedPymol.Source -cq -d "run view_6lzg_titan_vs_apbs.py; viewport $Width,$Height; png $rawName,$Width,$Height,300,0; quit"
    if ($LASTEXITCODE -ne 0) {
        throw "PyMOL render failed with exit code $LASTEXITCODE"
    }
} finally {
    Pop-Location
}
if (-not (Test-Path -LiteralPath $rawPath)) {
    throw "PyMOL did not write $rawPath"
}

$summary = Import-Csv -LiteralPath $summaryPath | Select-Object -First 1
$titanLimit = [double]$summary.titan_color_limit_kT_e
$apbsLimit = [double]$summary.apbs_color_limit_kT_e
if ([double]::IsNaN($titanLimit) -or [double]::IsInfinity($titanLimit) -or $titanLimit -le 0.0 -or
    [double]::IsNaN($apbsLimit) -or [double]::IsInfinity($apbsLimit) -or $apbsLimit -le 0.0) {
    throw 'Invalid color limits in the comparison summary.'
}

Add-Type -AssemblyName System.Drawing
Add-Type -ReferencedAssemblies System.Drawing -TypeDefinition @'
using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;

public static class TitanSurfaceFigureRenderer
{
    private static Rectangle FindContent(Bitmap image)
    {
        int minX = image.Width, minY = image.Height, maxX = -1, maxY = -1;
        for (int y = 0; y < image.Height; ++y)
        {
            for (int x = 0; x < image.Width; ++x)
            {
                Color c = image.GetPixel(x, y);
                if (c.A > 16 && (c.R < 246 || c.G < 246 || c.B < 246))
                {
                    if (x < minX) minX = x;
                    if (x > maxX) maxX = x;
                    if (y < minY) minY = y;
                    if (y > maxY) maxY = y;
                }
            }
        }
        if (maxX < minX || maxY < minY)
            throw new InvalidOperationException("Rendered image contains no visible molecular surface.");
        const int margin = 36;
        minX = Math.Max(0, minX - margin);
        minY = Math.Max(0, minY - margin);
        maxX = Math.Min(image.Width - 1, maxX + margin);
        maxY = Math.Min(image.Height - 1, maxY + margin);
        return Rectangle.FromLTRB(minX, minY, maxX + 1, maxY + 1);
    }

    private static void DrawCentered(Graphics g, string text, Font font, Brush brush,
                                     float centerX, float y)
    {
        SizeF size = g.MeasureString(text, font);
        g.DrawString(text, font, brush, centerX - 0.5f * size.Width, y);
    }

    private static void DrawColorBar(Graphics g, float x, float y, float width,
                                     float height, double limit, string label,
                                     Font labelFont, Font tickFont)
    {
        for (int i = 0; i < (int)width; ++i)
        {
            double t = 2.0 * i / Math.Max(1.0, width - 1.0) - 1.0;
            int r, green, b;
            if (t < 0.0)
            {
                r = 255;
                green = b = (int)Math.Round(255.0 * (1.0 + t));
            }
            else
            {
                r = green = (int)Math.Round(255.0 * (1.0 - t));
                b = 255;
            }
            using (var pen = new Pen(Color.FromArgb(r, green, b)))
                g.DrawLine(pen, x + i, y, x + i, y + height);
        }
        g.DrawRectangle(Pens.Black, x, y, width, height);
        DrawCentered(g, label, labelFont, Brushes.Black, x + 0.5f * width, y - 78.0f);
        string lo = String.Format(System.Globalization.CultureInfo.InvariantCulture, "-{0:0.##}", limit);
        string hi = String.Format(System.Globalization.CultureInfo.InvariantCulture, "+{0:0.##}", limit);
        g.DrawString(lo, tickFont, Brushes.Black, x, y + height + 8.0f);
        DrawCentered(g, "0", tickFont, Brushes.Black, x + 0.5f * width, y + height + 8.0f);
        SizeF hiSize = g.MeasureString(hi, tickFont);
        g.DrawString(hi, tickFont, Brushes.Black, x + width - hiSize.Width, y + height + 8.0f);
    }

    private static void DrawArrow(Graphics g, Pen pen, PointF from, PointF to)
    {
        g.DrawLine(pen, from, to);
        double angle = Math.Atan2(to.Y - from.Y, to.X - from.X);
        const float head = 18.0f;
        PointF a = new PointF(
            to.X - head * (float)Math.Cos(angle - Math.PI / 6.0),
            to.Y - head * (float)Math.Sin(angle - Math.PI / 6.0));
        PointF b = new PointF(
            to.X - head * (float)Math.Cos(angle + Math.PI / 6.0),
            to.Y - head * (float)Math.Sin(angle + Math.PI / 6.0));
        using (var brush = new SolidBrush(pen.Color))
            g.FillPolygon(brush, new PointF[] { to, a, b });
    }

    private static void DrawCharge(Graphics g, float x, float y, float radius,
                                   Color color, string sign, Font font)
    {
        using (var brush = new SolidBrush(color))
        {
            g.FillEllipse(brush, x - radius, y - radius, 2.0f * radius, 2.0f * radius);
            g.DrawEllipse(Pens.DimGray, x - radius, y - radius, 2.0f * radius, 2.0f * radius);
        }
        SizeF size = g.MeasureString(sign, font);
        g.DrawString(sign, font, Brushes.White, x - 0.5f * size.Width, y - 0.55f * size.Height);
    }

    private static void DrawCartoon(Graphics g, int width, Font panelFont,
                                    Font labelFont, Font smallFont,
                                    Font chargeFont)
    {
        const float y0 = 70.0f;
        g.DrawString("A", panelFont, Brushes.Black, 110.0f, y0);

        using (var linePen = new Pen(Color.FromArgb(70, 70, 70), 4.0f))
        using (var lightPen = new Pen(Color.FromArgb(135, 135, 135), 3.0f))
        using (var sesPen = new Pen(Color.FromArgb(30, 30, 30), 5.0f))
        using (var gridPen = new Pen(Color.FromArgb(165, 165, 165), 2.0f))
        using (var atomBrush = new SolidBrush(Color.FromArgb(220, 220, 220)))
        {
            float inputX = 430.0f;
            float cy = 370.0f;
            PointF[] atoms = {
                new PointF(inputX - 105, cy - 35), new PointF(inputX - 35, cy - 105),
                new PointF(inputX + 55, cy - 75), new PointF(inputX + 115, cy + 5),
                new PointF(inputX + 45, cy + 85), new PointF(inputX - 55, cy + 80),
                new PointF(inputX - 125, cy + 35), new PointF(inputX - 5, cy)
            };
            foreach (PointF p in atoms)
            {
                g.FillEllipse(atomBrush, p.X - 62, p.Y - 62, 124, 124);
                g.DrawEllipse(lightPen, p.X - 62, p.Y - 62, 124, 124);
            }
            DrawCharge(g, inputX - 80, cy - 30, 28, Color.FromArgb(230, 72, 72), "-", chargeFont);
            DrawCharge(g, inputX + 72, cy + 18, 28, Color.FromArgb(62, 92, 225), "+", chargeFont);
            DrawCentered(g, "same PQR charges", labelFont, Brushes.Black, inputX, 555.0f);

            float splitX = 920.0f;
            DrawArrow(g, linePen, new PointF(610, cy), new PointF(splitX, cy));
            DrawArrow(g, linePen, new PointF(splitX, cy), new PointF(1160, 235));
            DrawArrow(g, linePen, new PointF(splitX, cy), new PointF(1160, 500));

            RectangleF titanBox = new RectangleF(1180, 105, 690, 255);
            RectangleF apbsBox = new RectangleF(1180, 405, 690, 255);
            g.DrawRectangle(sesPen, titanBox.X, titanBox.Y, titanBox.Width, titanBox.Height);
            g.DrawRectangle(sesPen, apbsBox.X, apbsBox.Y, apbsBox.Width, apbsBox.Height);
            g.DrawString("TITAN", labelFont, Brushes.Black, 1215, 130);
            g.DrawString("analytic SES + node patches", smallFont, Brushes.Black, 1215, 190);
            g.DrawString("boundary evaluation", smallFont, Brushes.Black, 1215, 250);
            g.DrawString("APBS", labelFont, Brushes.Black, 1215, 430);
            g.DrawString("finite-difference grid", smallFont, Brushes.Black, 1215, 490);
            g.DrawString("surface sampling", smallFont, Brushes.Black, 1215, 550);

            float gridX = 1660.0f;
            float gridY = 470.0f;
            for (int i = 0; i <= 4; ++i)
            {
                g.DrawLine(gridPen, gridX + i * 34, gridY, gridX + i * 34, gridY + 136);
                g.DrawLine(gridPen, gridX, gridY + i * 34, gridX + 136, gridY + i * 34);
            }

            float outputX = 2420.0f;
            DrawArrow(g, linePen, new PointF(1870, 235), new PointF(outputX - 250, 235));
            DrawArrow(g, linePen, new PointF(1870, 535), new PointF(outputX - 250, 535));

            PointF[] upperSurface = {
                new PointF(outputX - 185, 230), new PointF(outputX - 120, 165),
                new PointF(outputX - 25, 180), new PointF(outputX + 45, 135),
                new PointF(outputX + 145, 190), new PointF(outputX + 185, 275),
                new PointF(outputX + 120, 340), new PointF(outputX + 20, 325),
                new PointF(outputX - 70, 350), new PointF(outputX - 165, 300)
            };
            PointF[] lowerSurface = new PointF[upperSurface.Length];
            for (int i = 0; i < upperSurface.Length; ++i)
                lowerSurface[i] = new PointF(upperSurface[i].X, upperSurface[i].Y + 300.0f);
            using (var neutralBrush = new SolidBrush(Color.FromArgb(235, 235, 235)))
            {
                g.FillPolygon(neutralBrush, upperSurface);
                g.FillPolygon(neutralBrush, lowerSurface);
            }
            g.DrawPolygon(sesPen, upperSurface);
            g.DrawPolygon(sesPen, lowerSurface);
            DrawCharge(g, outputX - 80, 245, 27, Color.FromArgb(230, 72, 72), "-", chargeFont);
            DrawCharge(g, outputX + 80, 255, 27, Color.FromArgb(62, 92, 225), "+", chargeFont);
            DrawCharge(g, outputX - 50, 545, 27, Color.FromArgb(230, 72, 72), "-", chargeFont);
            DrawCharge(g, outputX + 65, 530, 27, Color.FromArgb(62, 92, 225), "+", chargeFont);
            g.DrawString("boundary field", smallFont, Brushes.Black, outputX + 230, 245);
            g.DrawString("sampled field", smallFont, Brushes.Black, outputX + 230, 545);
        }
    }

    public static void Render(string rawPath, string outputPath,
                              double titanLimit, double apbsLimit)
    {
        using (var source = new Bitmap(rawPath))
        {
            Rectangle crop = FindContent(source);
            const int outputWidth = 3200;
            const int surfaceTop = 790;
            const int footer = 390;
            int surfaceHeight = (int)Math.Round(crop.Height * (outputWidth - 180.0) / crop.Width);
            int outputHeight = surfaceTop + surfaceHeight + footer;
            using (var output = new Bitmap(outputWidth, outputHeight, PixelFormat.Format32bppArgb))
            {
                output.SetResolution(600.0f, 600.0f);
                using (Graphics g = Graphics.FromImage(output))
                using (var panelFont = new Font("Arial", 72.0f, FontStyle.Bold, GraphicsUnit.Pixel))
                using (var labelFont = new Font("Arial", 58.0f, FontStyle.Regular, GraphicsUnit.Pixel))
                using (var tickFont = new Font("Arial", 48.0f, FontStyle.Regular, GraphicsUnit.Pixel))
                using (var cartoonLabelFont = new Font("Arial", 62.0f, FontStyle.Regular, GraphicsUnit.Pixel))
                using (var cartoonFont = new Font("Arial", 50.0f, FontStyle.Regular, GraphicsUnit.Pixel))
                using (var chargeFont = new Font("Arial", 34.0f, FontStyle.Bold, GraphicsUnit.Pixel))
                {
                    g.Clear(Color.White);
                    g.SmoothingMode = SmoothingMode.HighQuality;
                    g.InterpolationMode = InterpolationMode.HighQualityBicubic;
                    g.PixelOffsetMode = PixelOffsetMode.HighQuality;
                    DrawCartoon(g, outputWidth, panelFont, cartoonLabelFont,
                                cartoonFont, chargeFont);

                    Rectangle destination = new Rectangle(90, surfaceTop,
                                                          outputWidth - 180, surfaceHeight);
                    g.DrawImage(source, destination, crop, GraphicsUnit.Pixel);

                    float leftCenter = 0.29f * outputWidth;
                    float rightCenter = 0.71f * outputWidth;
                    g.DrawString("B", panelFont, Brushes.Black, 110.0f, surfaceTop - 88.0f);
                    g.DrawString("C", panelFont, Brushes.Black, 0.51f * outputWidth, surfaceTop - 88.0f);

                    float barY = surfaceTop + surfaceHeight + 160.0f;
                    DrawColorBar(g, leftCenter - 510.0f, barY, 1020.0f, 46.0f,
                                 titanLimit, "TITAN boundary potential (kT/e)",
                                 labelFont, tickFont);
                    DrawColorBar(g, rightCenter - 510.0f, barY, 1020.0f, 46.0f,
                                 apbsLimit, "APBS grid potential sampled on SES (kT/e)",
                                 labelFont, tickFont);
                }
                output.Save(outputPath, ImageFormat.Png);
            }
        }
    }
}
'@

[TitanSurfaceFigureRenderer]::Render($rawPath, $OutputPath, $titanLimit, $apbsLimit)
Write-Host "Publication figure: $OutputPath"
Write-Host ("Independent robust color limits: TITAN +/-{0:G5} kT/e; APBS +/-{1:G5} kT/e" -f $titanLimit, $apbsLimit)
